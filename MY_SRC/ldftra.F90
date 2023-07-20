MODULE ldftra
   !!======================================================================
   !!                       ***  MODULE  ldftra  ***
   !! Ocean physics:  lateral diffusivity coefficients
   !!=====================================================================
   !! History :       ! 1997-07  (G. Madec)  from inimix.F split in 2 routines
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2005-11  (G. Madec)
   !!            3.7  ! 2013-12  (F. Lemarie, G. Madec)  restructuration/simplification of aht/aeiv specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_tra_init : initialization, namelist read, and parameters control
   !!   ldf_tra      : update lateral eddy diffusivity coefficients at each time step
   !!   ldf_eiv_init : initialization of the eiv coeff. from namelist choices
   !!   ldf_eiv      : time evolution of the eiv coefficients (function of the growth rate of baroclinic instability)
   !!   ldf_eiv_trp  : add to the input ocean transport the contribution of the EIV parametrization
   !!   ldf_eiv_dia  : diagnose the eddy induced velocity from the eiv streamfunction
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldfslp          ! lateral diffusion: slope of iso-neutral surfaces
   USE ldfc1d_c2d      ! lateral diffusion: 1D & 2D cases
   USE diaptr
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module for ehanced bottom friction file
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_tra_init   ! called by nemogcm.F90
   PUBLIC   ldf_tra        ! called by step.F90
   PUBLIC   ldf_eiv_init   ! called by nemogcm.F90
   PUBLIC   ldf_eiv        ! called by step.F90
   PUBLIC   ldf_eiv_trp    ! called by traadv.F90
   PUBLIC   ldf_eiv_dia    ! called by traldf_iso and traldf_iso_triad.F90

   !                                   !!* Namelist namtra_ldf : lateral mixing on tracers *
   !                                    != Operator type =!
   LOGICAL , PUBLIC ::   ln_traldf_OFF       !: no operator: No explicit diffusion
   LOGICAL , PUBLIC ::   ln_traldf_lap       !: laplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_blp       !: bilaplacian operator
   !                                    != Direction of action =!
   LOGICAL , PUBLIC ::   ln_traldf_lev       !: iso-level direction
   LOGICAL , PUBLIC ::   ln_traldf_hor       !: horizontal (geopotential) direction
!  LOGICAL , PUBLIC ::   ln_traldf_iso       !: iso-neutral direction                    (see ldfslp)
   !		       	                      != iso-neutral options =!
!  LOGICAL , PUBLIC ::   ln_traldf_triad     !: griffies triad scheme                    (see ldfslp)
   LOGICAL , PUBLIC ::   ln_traldf_msc       !: Method of Stabilizing Correction
!  LOGICAL , PUBLIC ::   ln_triad_iso        !: pure horizontal mixing in ML             (see ldfslp)
!  LOGICAL , PUBLIC ::   ln_botmix_triad     !: mixing on bottom                         (see ldfslp)
!  REAL(wp), PUBLIC ::   rn_sw_triad         !: =1/0 switching triad / all 4 triads used (see ldfslp)
!  REAL(wp), PUBLIC ::   rn_slpmax           !: slope limit                              (see ldfslp)
   !                                    !=  Coefficients =!
   INTEGER , PUBLIC ::   nn_aht_ijk_t        !: choice of time & space variations of the lateral eddy diffusivity coef.
   !                                            !  time invariant coefficients:  aht_0 = 1/2  Ud*Ld   (lap case)
   !                                            !                                bht_0 = 1/12 Ud*Ld^3 (blp case)
   REAL(wp), PUBLIC ::      rn_Ud               !: lateral diffusive velocity  [m/s]
   REAL(wp), PUBLIC ::      rn_Ld               !: lateral diffusive length    [m]

   !                                   !!* Namelist namtra_eiv : eddy induced velocity param. *
   !                                    != Use/diagnose eiv =!
   !LOGICAL , PUBLIC ::   ln_ldfeiv           !: eddy induced velocity flag
   LOGICAL , PUBLIC ::   ln_ldfeiv_dia       !: diagnose & output eiv streamfunction and velocity (IOM)
   !                                    != Coefficients =!
   INTEGER , PUBLIC ::   nn_aei_ijk_t        !: choice of time/space variation of the eiv coeff.
   REAL(wp), PUBLIC ::      rn_Ue               !: lateral diffusive velocity  [m/s]
   REAL(wp), PUBLIC ::      rn_Le               !: lateral diffusive length    [m]

   !                                  ! Flag to control the type of lateral diffusive operator
   INTEGER, PARAMETER, PUBLIC ::   np_ERROR  =-10   ! error in specification of lateral diffusion
   INTEGER, PARAMETER, PUBLIC ::   np_no_ldf = 00   ! without operator (i.e. no lateral diffusive trend)
   !                          !!      laplacian     !    bilaplacian    !
   INTEGER, PARAMETER, PUBLIC ::   np_lap    = 10   ,   np_blp    = 20  ! iso-level operator
   INTEGER, PARAMETER, PUBLIC ::   np_lap_i  = 11   ,   np_blp_i  = 21  ! standard iso-neutral or geopotential operator
   INTEGER, PARAMETER, PUBLIC ::   np_lap_it = 12   ,   np_blp_it = 22  ! triad    iso-neutral or geopotential operator

   INTEGER , PUBLIC ::   nldf_tra      = 0         !: type of lateral diffusion used defined from ln_traldf_... (namlist logicals)
   LOGICAL , PUBLIC ::   l_ldftra_time = .FALSE.   !: flag for time variation of the lateral eddy diffusivity coef.
   LOGICAL , PUBLIC ::   l_ldfeiv_time = .FALSE.   !: flag for time variation of the eiv coef.

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahtu, ahtv   !: eddy diffusivity coef. at U- and V-points   [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aeiu, aeiv   !: eddy induced velocity coeff.                [m2/s]

   REAL(wp) ::   aht0, aei0               ! constant eddy coefficients (deduced from namelist values)                     [m2/s]
   REAL(wp) ::   r1_2  = 0.5_wp           ! =1/2
   REAL(wp) ::   r1_4  = 0.25_wp          ! =1/4
   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! =1/12

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: ldftra.F90 15475 2021-11-05 14:14:45Z cdllod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_tra_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_init  ***
      !!
      !! ** Purpose :   initializations of the tracer lateral mixing coeff.
      !!
      !! ** Method  : * the eddy diffusivity coef. specification depends on:
      !!
      !!    ln_traldf_lap = T     laplacian operator
      !!    ln_traldf_blp = T   bilaplacian operator
      !!
      !!    nn_aht_ijk_t  =  0 => = constant
      !!                  !
      !!                  = 10 => = F(z) : constant with a reduction of 1/4 with depth
      !!                  !
      !!                  =-20 => = F(i,j)   = shape read in 'eddy_diffusivity_2D.nc' file
      !!                  = 20    = F(i,j)   = F(e1,e2) or F(e1^3,e2^3) (lap or bilap case)
      !!                  = 21    = F(i,j,t) = F(growth rate of baroclinic instability)
      !!                  !
      !!                  =-30 => = F(i,j,k)   = shape read in 'eddy_diffusivity_3D.nc' file
      !!                  = 30    = F(i,j,k)   = 2D (case 20) + decrease with depth (case 10)
      !!                  = 31    = F(i,j,k,t) = F(local velocity) (  1/2  |u|e     laplacian operator
      !!                                                           or 1/12 |u|e^3 bilaplacian operator )
      !!              * initialisation of the eddy induced velocity coefficient by a call to ldf_eiv_init
      !!
      !! ** action  : ahtu, ahtv initialized one for all or l_ldftra_time set to true
      !!              aeiu, aeiv initialized one for all or l_ldfeiv_time set to true
      !!----------------------------------------------------------------------
      INTEGER  ::   jk                             ! dummy loop indices
      INTEGER  ::   ioptio, ierr, inum, ios, inn   ! local integer
      REAL(wp) ::   zah_max, zUfac                 !   -      -
      CHARACTER(len=5) ::   cl_Units               ! units (m2/s or m4/s)
      !!
      NAMELIST/namtra_ldf/ ln_traldf_OFF, ln_traldf_lap  , ln_traldf_blp  ,   &   ! type of operator
         &                 ln_traldf_lev, ln_traldf_hor  , ln_traldf_triad,   &   ! acting direction of the operator
         &                 ln_traldf_iso, ln_traldf_msc  ,  rn_slpmax     ,   &   ! option for iso-neutral operator
         &                 ln_triad_iso , ln_botmix_triad, rn_sw_triad    ,   &   ! option for triad operator
         &                 nn_aht_ijk_t , rn_Ud          , rn_Ld                  ! lateral eddy coefficient
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra_init : lateral tracer diffusion'
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF

      !
      !  Choice of lateral tracer physics
      ! =================================
      !
      READ  ( numnam_ref, namtra_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_ldf in reference namelist' )
      READ  ( numnam_cfg, namtra_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_ldf in configuration namelist' )
      IF(lwm) WRITE( numond, namtra_ldf )
      !
      IF(lwp) THEN                      ! control print
         WRITE(numout,*) '   Namelist : namtra_ldf --- lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      type :'
         WRITE(numout,*) '         no explicit diffusion                   ln_traldf_OFF   = ', ln_traldf_OFF
         WRITE(numout,*) '         laplacian operator                      ln_traldf_lap   = ', ln_traldf_lap
         WRITE(numout,*) '         bilaplacian operator                    ln_traldf_blp   = ', ln_traldf_blp
         WRITE(numout,*) '      direction of action :'
         WRITE(numout,*) '         iso-level                               ln_traldf_lev   = ', ln_traldf_lev
         WRITE(numout,*) '         horizontal (geopotential)               ln_traldf_hor   = ', ln_traldf_hor
         WRITE(numout,*) '         iso-neutral Madec operator              ln_traldf_iso   = ', ln_traldf_iso
         WRITE(numout,*) '         iso-neutral triad operator              ln_traldf_triad = ', ln_traldf_triad
         WRITE(numout,*) '            use the Method of Stab. Correction   ln_traldf_msc   = ', ln_traldf_msc
         WRITE(numout,*) '            maximum isoppycnal slope             rn_slpmax       = ', rn_slpmax
         WRITE(numout,*) '            pure lateral mixing in ML            ln_triad_iso    = ', ln_triad_iso
         WRITE(numout,*) '            switching triad or not               rn_sw_triad     = ', rn_sw_triad
         WRITE(numout,*) '            lateral mixing on bottom             ln_botmix_triad = ', ln_botmix_triad
         WRITE(numout,*) '      coefficients :'
         WRITE(numout,*) '         type of time-space variation            nn_aht_ijk_t    = ', nn_aht_ijk_t
         WRITE(numout,*) '            lateral diffusive velocity (if cst)  rn_Ud           = ', rn_Ud, ' m/s'
         WRITE(numout,*) '            lateral diffusive length   (if cst)  rn_Ld           = ', rn_Ld, ' m'
      ENDIF
      !
      !
      ! Operator and its acting direction   (set nldf_tra)
      ! =================================
      !
      nldf_tra = np_ERROR
      ioptio   = 0
      IF( ln_traldf_OFF ) THEN   ;   nldf_tra = np_no_ldf   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_traldf_lap ) THEN   ;                              ioptio = ioptio + 1   ;   ENDIF
      IF( ln_traldf_blp ) THEN   ;                              ioptio = ioptio + 1   ;   ENDIF
      IF( ioptio /=  1  )   CALL ctl_stop( 'tra_ldf_init: use ONE of the 3 operator options (NONE/lap/blp)' )
      !
      IF( .NOT.ln_traldf_OFF ) THEN    !==  direction ==>> type of operator  ==!
         ioptio = 0
         IF( ln_traldf_lev   )   ioptio = ioptio + 1
         IF( ln_traldf_hor   )   ioptio = ioptio + 1
         IF( ln_traldf_iso   )   ioptio = ioptio + 1
         IF( ln_traldf_triad )   ioptio = ioptio + 1
         IF( ioptio /=  1  )   CALL ctl_stop( 'tra_ldf_init: use ONE direction (level/hor/iso/triad)' )
         !
         !                                ! defined the type of lateral diffusion from ln_traldf_... logicals
         ierr = 0
         IF ( ln_traldf_lap ) THEN        ! laplacian operator
            IF ( ln_zco ) THEN                  ! z-coordinate
               IF ( ln_traldf_lev   )   nldf_tra = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_hor   )   nldf_tra = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_iso   )   nldf_tra = np_lap_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf_tra = np_lap_it  ! iso-neutral: triad     (   rotation)
            ENDIF
            IF ( ln_zps ) THEN                  ! z-coordinate with partial step
               IF ( ln_traldf_lev   )   ierr     = 1          ! iso-level not allowed
               IF ( ln_traldf_hor   )   nldf_tra = np_lap     ! horizontal             (no rotation)
               IF ( ln_traldf_iso   )   nldf_tra = np_lap_i   ! iso-neutral: standard     (rotation)
               IF ( ln_traldf_triad )   nldf_tra = np_lap_it  ! iso-neutral: triad        (rotation)
            ENDIF
            IF ( ln_sco ) THEN                  ! s-coordinate
               IF ( ln_traldf_lev   )   nldf_tra = np_lap     ! iso-level              (no rotation)
               IF ( ln_traldf_hor   )   nldf_tra = np_lap_i   ! horizontal             (   rotation)
               IF ( ln_traldf_iso   )   nldf_tra = np_lap_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf_tra = np_lap_it  ! iso-neutral: triad     (   rotation)
            ENDIF
         ENDIF
         !
         IF( ln_traldf_blp ) THEN         ! bilaplacian operator
            IF ( ln_zco ) THEN                  ! z-coordinate
               IF ( ln_traldf_lev   )   nldf_tra = np_blp     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_hor   )   nldf_tra = np_blp     ! iso-level = horizontal (no rotation)
               IF ( ln_traldf_iso   )   nldf_tra = np_blp_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf_tra = np_blp_it  ! iso-neutral: triad     (   rotation)
            ENDIF
            IF ( ln_zps ) THEN                  ! z-coordinate with partial step
               IF ( ln_traldf_lev   )   ierr     = 1          ! iso-level not allowed
               IF ( ln_traldf_hor   )   nldf_tra = np_blp     ! horizontal             (no rotation)
               IF ( ln_traldf_iso   )   nldf_tra = np_blp_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf_tra = np_blp_it  ! iso-neutral: triad     (   rotation)
            ENDIF
            IF ( ln_sco ) THEN                  ! s-coordinate
               IF ( ln_traldf_lev   )   nldf_tra = np_blp     ! iso-level              (no rotation)
               IF ( ln_traldf_hor   )   nldf_tra = np_blp_it  ! horizontal             (   rotation)
               IF ( ln_traldf_iso   )   nldf_tra = np_blp_i   ! iso-neutral: standard  (   rotation)
               IF ( ln_traldf_triad )   nldf_tra = np_blp_it  ! iso-neutral: triad     (   rotation)
            ENDIF
         ENDIF
         IF ( ierr == 1 )   CALL ctl_stop( 'iso-level in z-partial step, not allowed' )
      ENDIF
      !
      IF( ln_isfcav .AND. ln_traldf_triad )   CALL ctl_stop( ' ice shelf cavity and traldf_triad not tested' )
           !
      IF(  nldf_tra == np_lap_i .OR. nldf_tra == np_lap_it .OR. &
         & nldf_tra == np_blp_i .OR. nldf_tra == np_blp_it  )   l_ldfslp = .TRUE.    ! slope of neutral surfaces required
      !
      IF( ln_traldf_blp .AND. ( ln_traldf_iso .OR. ln_traldf_triad) ) THEN     ! iso-neutral bilaplacian need MSC
         IF( .NOT.ln_traldf_msc )   CALL ctl_stop( 'tra_ldf_init: iso-neutral bilaplacian requires ln_traldf_msc=.true.' )
      ENDIF
      !
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nldf_tra )
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral diffusion'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   laplacian iso-level operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   Rotated laplacian operator (standard)'
         CASE( np_lap_it )   ;   WRITE(numout,*) '   ==>>>   Rotated laplacian operator (triad)'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   bilaplacian iso-level operator'
         CASE( np_blp_i  )   ;   WRITE(numout,*) '   ==>>>   Rotated bilaplacian operator (standard)'
         CASE( np_blp_it )   ;   WRITE(numout,*) '   ==>>>   Rotated bilaplacian operator (triad)'
         END SELECT
         WRITE(numout,*)
      ENDIF

      !
      !  Space/time variation of eddy coefficients
      ! ===========================================
      !
      l_ldftra_time = .FALSE.                ! no time variation except in case defined below
      !
      IF( ln_traldf_OFF ) THEN               !== no explicit diffusive operator  ==!
         !
         IF(lwp) WRITE(numout,*) '   ==>>>   No diffusive operator selected. ahtu and ahtv are not allocated'
         RETURN
         !
      ELSE                                   !==  a lateral diffusion operator is used  ==!
         !
         !                                         ! allocate the aht arrays
         ALLOCATE( ahtu(jpi,jpj,jpk) , ahtv(jpi,jpj,jpk) , STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_tra_init: failed to allocate arrays')
         !
         ahtu(:,:,jpk) = 0._wp                     ! last level always 0
         ahtv(:,:,jpk) = 0._wp
         !.
         !                                         ! value of lap/blp eddy mixing coef.
         IF(     ln_traldf_lap ) THEN   ;   zUfac = r1_2 *rn_Ud   ;   inn = 1   ;   cl_Units = ' m2/s'   !   laplacian
         ELSEIF( ln_traldf_blp ) THEN   ;   zUfac = r1_12*rn_Ud   ;   inn = 3   ;   cl_Units = ' m4/s'   ! bilaplacian
         ENDIF
         aht0    = zUfac *    rn_Ld**inn              ! mixing coefficient
         zah_max = zUfac * (ra*rad)**inn              ! maximum reachable coefficient (value at the Equator for e1=1 degree)
         !
         !
         SELECT CASE(  nn_aht_ijk_t  )             !* Specification of space-time variations of ahtu, ahtv
         !
         CASE(   0  )      !==  constant  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = constant = ', aht0, cl_Units
            ahtu(:,:,1:jpkm1) = aht0
            ahtv(:,:,1:jpkm1) = aht0
            !
         CASE(  10  )      !==  fixed profile  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F( depth )'
            IF(lwp) WRITE(numout,*) '           surface eddy diffusivity = constant = ', aht0, cl_Units
            ahtu(:,:,1) = aht0                        ! constant surface value
            ahtv(:,:,1) = aht0
            CALL ldf_c1d( 'TRA', ahtu(:,:,1), ahtv(:,:,1), ahtu, ahtv )
            !
         CASE ( -20 )      !== fixed horizontal shape and magnitude read in file  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F(i,j) read in eddy_diffusivity.nc file'
            CALL iom_open( 'eddy_diffusivity_2D.nc', inum )
            CALL iom_get ( inum, jpdom_global, 'ahtu_2D', ahtu(:,:,1), cd_type = 'U', psgn = 1._wp )
            CALL iom_get ( inum, jpdom_global, 'ahtv_2D', ahtv(:,:,1), cd_type = 'V', psgn = 1._wp )
            CALL iom_close( inum )
            DO jk = 2, jpkm1
               ahtu(:,:,jk) = ahtu(:,:,1)
               ahtv(:,:,jk) = ahtv(:,:,1)
            END DO
            !
         CASE(  20  )      !== fixed horizontal shape  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F( e1, e2 ) or F( e1^3, e2^3 ) (lap or blp case)'
            IF(lwp) WRITE(numout,*) '           using a fixed diffusive velocity = ', rn_Ud,' m/s   and   Ld = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
            CALL ldf_c2d( 'TRA', zUfac      , inn        , ahtu, ahtv )    ! value proportional to scale factor^inn
            !
         CASE(  21  )      !==  time varying 2D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F( latitude, longitude, time )'
            IF(lwp) WRITE(numout,*) '                            = F( growth rate of baroclinic instability )'
            IF(lwp) WRITE(numout,*) '            min value = 0.2 * aht0 (with aht0= 1/2 rn_Ud*rn_Ld)'
            IF(lwp) WRITE(numout,*) '            max value =       aei0 (with aei0=1/2 rn_Ue*Le  increased to aht0 within 20N-20S'
            !
            l_ldftra_time = .TRUE.     ! will be calculated by call to ldf_tra routine in step.F90
            !
            IF( ln_traldf_blp )   CALL ctl_stop( 'ldf_tra_init: aht=F( growth rate of baroc. insta .)',   &
               &                                 '              incompatible with bilaplacian operator' )
            !
         CASE( -30  )      !== fixed 3D shape read in file  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F(i,j,k) read in eddy_diffusivity.nc file'
            CALL iom_open( 'eddy_diffusivity_3D.nc', inum )
            CALL iom_get ( inum, jpdom_global, 'ahtu_3D', ahtu, cd_type = 'U', psgn = 1._wp )
            CALL iom_get ( inum, jpdom_global, 'ahtv_3D', ahtv, cd_type = 'V', psgn = 1._wp )
            CALL iom_close( inum )
            !
         CASE(  30  )      !==  fixed 3D shape  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F( latitude, longitude, depth )'
            IF(lwp) WRITE(numout,*) '           using a fixed diffusive velocity = ', rn_Ud,' m/s   and   Ld = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
            CALL ldf_c2d( 'TRA', zUfac      , inn        , ahtu, ahtv )    ! surface value proportional to scale factor^inn
            CALL ldf_c1d( 'TRA', ahtu(:,:,1), ahtv(:,:,1), ahtu, ahtv )    ! reduction with depth
            !
         CASE(  31  )      !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy diffusivity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the velocity : 1/2 |u|e or 1/12 |u|e^3'
            !
            l_ldftra_time = .TRUE.     ! will be calculated by call to ldf_tra routine in step.F90
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_tra_init: wrong choice for nn_aht_ijk_t, the type of space-time variation of aht')
         END SELECT
         !
         IF( .NOT.l_ldftra_time ) THEN             !* No time variation
            IF(     ln_traldf_lap ) THEN                 !   laplacian operator (mask only)
               ahtu(:,:,1:jpkm1) =       ahtu(:,:,1:jpkm1)   * umask(:,:,1:jpkm1)
               ahtv(:,:,1:jpkm1) =       ahtv(:,:,1:jpkm1)   * vmask(:,:,1:jpkm1)
            ELSEIF( ln_traldf_blp ) THEN                 ! bilaplacian operator (square root + mask)
               ahtu(:,:,1:jpkm1) = SQRT( ahtu(:,:,1:jpkm1) ) * umask(:,:,1:jpkm1)
               ahtv(:,:,1:jpkm1) = SQRT( ahtv(:,:,1:jpkm1) ) * vmask(:,:,1:jpkm1)
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_tra_init


   SUBROUTINE ldf_tra( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra  ***
      !!
      !! ** Purpose :   update at kt the tracer lateral mixing coeff. (aht and aeiv)
      !!
      !! ** Method  : * time varying eddy diffusivity coefficients:
      !!
      !!    nn_aei_ijk_t = 21    aeiu, aeiv = F(i,j,  t) = F(growth rate of baroclinic instability)
      !!                                                   with a reduction to 0 in vicinity of the Equator
      !!    nn_aht_ijk_t = 21    ahtu, ahtv = F(i,j,  t) = F(growth rate of baroclinic instability)
      !!
      !!                 = 31    ahtu, ahtv = F(i,j,k,t) = F(local velocity) (  |u|e  /12   laplacian operator
      !!                                                                     or |u|e^3/12 bilaplacian operator )
      !!
      !!              * time varying EIV coefficients: call to ldf_eiv routine
      !!
      !! ** action  :   ahtu, ahtv   update at each time step
      !!                aeiu, aeiv      -       -     -    -   (if ln_ldfeiv=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zaht, zahf, zaht_min, zDaht, z1_f20   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( ln_ldfeiv .AND. nn_aei_ijk_t == 21 ) THEN       ! eddy induced velocity coefficients
         !                                ! =F(growth rate of baroclinic instability)
         !                                ! max value aeiv_0 ; decreased to 0 within 20N-20S
         CALL ldf_eiv( kt, aei0, aeiu, aeiv, Kmm )
      ENDIF
      !
      SELECT CASE(  nn_aht_ijk_t  )       ! Eddy diffusivity coefficients
      !
      CASE(  21  )       !==  time varying 2D field  ==!   = F( growth rate of baroclinic instability )
         !                                             !   min value 0.2*aht0
         !                                             !   max value aht0 (aei0 if nn_aei_ijk_t=21)
         !                                             !   increase to aht0 within 20N-20S
         IF( ln_ldfeiv .AND. nn_aei_ijk_t == 21 ) THEN   ! use the already computed aei.
            ahtu(:,:,1) = aeiu(:,:,1)
            ahtv(:,:,1) = aeiv(:,:,1)
         ELSE                                            ! compute aht.
            CALL ldf_eiv( kt, aht0, ahtu, ahtv, Kmm )
         ENDIF
         !
         z1_f20   = 1._wp / (  2._wp * omega * SIN( rad * 20._wp )  )   ! 1 / ff(20 degrees)
         zaht_min = 0.2_wp * aht0                                       ! minimum value for aht
         zDaht    = aht0 - zaht_min
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !!gm CAUTION : here we assume lat/lon grid in 20deg N/S band (like all ORCA cfg)
            !!     ==>>>   The Coriolis value is identical for t- & u_points, and for v- and f-points
            zaht = ( 1._wp -  MIN( 1._wp , ABS( ff_t(ji,jj) * z1_f20 ) ) ) * zDaht
            zahf = ( 1._wp -  MIN( 1._wp , ABS( ff_f(ji,jj) * z1_f20 ) ) ) * zDaht
            ahtu(ji,jj,1) = (  MAX( zaht_min, ahtu(ji,jj,1) ) + zaht  )     ! min value zaht_min
            ahtv(ji,jj,1) = (  MAX( zaht_min, ahtv(ji,jj,1) ) + zahf  )     ! increase within 20S-20N
         END_2D
         DO jk = 1, jpkm1                             ! deeper value = surface value + mask for all levels
            ahtu(:,:,jk) = ahtu(:,:,1) * umask(:,:,jk)
            ahtv(:,:,jk) = ahtv(:,:,1) * vmask(:,:,jk)
         END DO
         !
      CASE(  31  )       !==  time varying 3D field  ==!   = F( local velocity )
         IF( ln_traldf_lap     ) THEN          !   laplacian operator |u| e /12
            DO jk = 1, jpkm1
               ahtu(:,:,jk) = ABS( uu(:,:,jk,Kbb) ) * e1u(:,:) * r1_12   ! n.b. uu,vv are masked
               ahtv(:,:,jk) = ABS( vv(:,:,jk,Kbb) ) * e2v(:,:) * r1_12
            END DO
         ELSEIF( ln_traldf_blp ) THEN      ! bilaplacian operator      sqrt( |u| e^3 /12 ) = sqrt( |u| e /12 ) * e
            DO jk = 1, jpkm1
               ahtu(:,:,jk) = SQRT(  ABS( uu(:,:,jk,Kbb) ) * e1u(:,:) * r1_12  ) * e1u(:,:)
               ahtv(:,:,jk) = SQRT(  ABS( vv(:,:,jk,Kbb) ) * e2v(:,:) * r1_12  ) * e2v(:,:)
            END DO
         ENDIF
         !
      END SELECT
      !
      CALL iom_put( "ahtu_2d", ahtu(:,:,1) )   ! surface u-eddy diffusivity coeff.
      CALL iom_put( "ahtv_2d", ahtv(:,:,1) )   ! surface v-eddy diffusivity coeff.
      CALL iom_put( "ahtu_3d", ahtu(:,:,:) )   ! 3D      u-eddy diffusivity coeff.
      CALL iom_put( "ahtv_3d", ahtv(:,:,:) )   ! 3D      v-eddy diffusivity coeff.
      !
      IF( ln_ldfeiv ) THEN
        CALL iom_put( "aeiu_2d", aeiu(:,:,1) )   ! surface u-EIV coeff.
        CALL iom_put( "aeiv_2d", aeiv(:,:,1) )   ! surface v-EIV coeff.
        CALL iom_put( "aeiu_3d", aeiu(:,:,:) )   ! 3D      u-EIV coeff.
        CALL iom_put( "aeiv_3d", aeiv(:,:,:) )   ! 3D      v-EIV coeff.
      ENDIF
      !
   END SUBROUTINE ldf_tra


   SUBROUTINE ldf_eiv_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv_init  ***
      !!
      !! ** Purpose :   initialization of the eiv coeff. from namelist choices.
      !!
      !! ** Method  :   the eiv diffusivity coef. specification depends on:
      !!    nn_aei_ijk_t  =  0 => = constant
      !!                  !
      !!                  = 10 => = F(z) : constant with a reduction of 1/4 with depth
      !!                  !
      !!                  =-20 => = F(i,j)   = shape read in 'eddy_induced_velocity_2D.nc' file
      !!                  = 20    = F(i,j)   = F(e1,e2) or F(e1^3,e2^3) (lap or bilap case)
      !!                  = 21    = F(i,j,t) = F(growth rate of baroclinic instability)
      !!                  !
      !!                  =-30 => = F(i,j,k)   = shape read in 'eddy_induced_velocity_3D.nc' file
      !!                  = 30    = F(i,j,k)   = 2D (case 20) + decrease with depth (case 10)
      !!
      !! ** Action  :   aeiu , aeiv   :  initialized one for all or l_ldftra_time set to true
      !!                l_ldfeiv_time : =T if EIV coefficients vary with time
      !!----------------------------------------------------------------------
      INTEGER  ::   jk                     ! dummy loop indices
      INTEGER  ::   ierr, inum, ios, inn   ! local integer
      REAL(wp) ::   zah_max, zUfac         !   -   scalar
      !!
      NAMELIST/namtra_eiv/ ln_ldfeiv   , ln_ldfeiv_dia,   &   ! eddy induced velocity (eiv)
         &                 nn_aei_ijk_t, rn_Ue, rn_Le         ! eiv  coefficient
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_eiv_init : eddy induced velocity parametrization'
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !
      READ  ( numnam_ref, namtra_eiv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_eiv in reference namelist' )
      !
      READ  ( numnam_cfg, namtra_eiv, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_eiv in configuration namelist' )
      IF(lwm)  WRITE ( numond, namtra_eiv )

      IF(lwp) THEN                      ! control print
         WRITE(numout,*) '   Namelist namtra_eiv : '
         WRITE(numout,*) '      Eddy Induced Velocity (eiv) param.         ln_ldfeiv     = ', ln_ldfeiv
         WRITE(numout,*) '      eiv streamfunction & velocity diag.        ln_ldfeiv_dia = ', ln_ldfeiv_dia
         WRITE(numout,*) '      coefficients :'
         WRITE(numout,*) '         type of time-space variation            nn_aei_ijk_t  = ', nn_aei_ijk_t
         WRITE(numout,*) '         lateral diffusive velocity (if cst)     rn_Ue         = ', rn_Ue, ' m/s'
         WRITE(numout,*) '         lateral diffusive length   (if cst)     rn_Le         = ', rn_Le, ' m'
         WRITE(numout,*)
      ENDIF
      !
      l_ldfeiv_time = .FALSE.       ! no time variation except in case defined below
      !
      !
      IF( .NOT.ln_ldfeiv ) THEN     !== Parametrization not used  ==!
         !
         IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity param is NOT used'
         ln_ldfeiv_dia = .FALSE.
         !
      ELSE                          !== use the parametrization  ==!
         !
         IF(lwp) WRITE(numout,*) '   ==>>>   use eddy induced velocity parametrization'
         IF(lwp) WRITE(numout,*)
         !
         IF( ln_traldf_blp )   CALL ctl_stop( 'ldf_eiv_init: eddy induced velocity ONLY with laplacian diffusivity' )
         !
         IF( .NOT.( ln_traldf_iso .OR. ln_traldf_triad ) )   &
           &                  CALL ctl_stop( 'ln_ldfeiv=T requires iso-neutral laplacian diffusion' )
         !                                != allocate the aei arrays
         ALLOCATE( aeiu(jpi,jpj,jpk), aeiv(jpi,jpj,jpk), STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ldf_eiv: failed to allocate arrays')
         !
         !                                != Specification of space-time variations of eaiu, aeiv
         !
         aeiu(:,:,jpk) = 0._wp               ! last level always 0
         aeiv(:,:,jpk) = 0._wp
         !                                   ! value of EIV coef. (laplacian operator)
         zUfac = r1_2 *rn_Ue                    ! velocity factor
         inn = 1                                ! L-exponent
         aei0    = zUfac *    rn_Le**inn        ! mixing coefficient
         zah_max = zUfac * (ra*rad)**inn        ! maximum reachable coefficient (value at the Equator)

         SELECT CASE( nn_aei_ijk_t )         !* Specification of space-time variations
         !
         CASE(   0  )                        !--  constant  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = constant = ', aei0, ' m2/s'
            aeiu(:,:,1:jpkm1) = aei0
            aeiv(:,:,1:jpkm1) = aei0
            !
         CASE(  10  )                        !--  fixed profile  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = F( depth )'
            IF(lwp) WRITE(numout,*) '           surface eddy diffusivity = constant = ', aht0, ' m2/s'
            aeiu(:,:,1) = aei0                  ! constant surface value
            aeiv(:,:,1) = aei0
            CALL ldf_c1d( 'TRA', aeiu(:,:,1), aeiv(:,:,1), aeiu, aeiv )
            !
         CASE ( -20 )                        !--  fixed horizontal shape read in file  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = F(i,j) read in eddy_diffusivity_2D.nc file'
            CALL iom_open ( 'eddy_induced_velocity_2D.nc', inum )
            CALL iom_get  ( inum, jpdom_global, 'aeiu', aeiu(:,:,1), cd_type = 'U', psgn = 1._wp )
            CALL iom_get  ( inum, jpdom_global, 'aeiv', aeiv(:,:,1), cd_type = 'V', psgn = 1._wp )
            CALL iom_close( inum )
            DO jk = 2, jpkm1
               aeiu(:,:,jk) = aeiu(:,:,1)
               aeiv(:,:,jk) = aeiv(:,:,1)
            END DO
            !
         CASE(  20  )                        !--  fixed horizontal shape  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = F( e1, e2 )'
            IF(lwp) WRITE(numout,*) '           using a fixed diffusive velocity = ', rn_Ue, ' m/s   and   Le = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, ' m2/s   for e1=1°)'
            CALL ldf_c2d( 'TRA', zUfac      , inn        , aeiu, aeiv )    ! value proportional to scale factor^inn
            !
         CASE(  21  )                        !--  time varying 2D field  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = F( latitude, longitude, time )'
            IF(lwp) WRITE(numout,*) '                                       = F( growth rate of baroclinic instability )'
            IF(lwp) WRITE(numout,*) '           maximum allowed value: aei0 = ', aei0, ' m2/s'
            !
            l_ldfeiv_time = .TRUE.     ! will be calculated by call to ldf_tra routine in step.F90
            !
         CASE( -30  )                        !-- fixed 3D shape read in file  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = F(i,j,k) read in eddy_diffusivity_3D.nc file'
            CALL iom_open ( 'eddy_induced_velocity_3D.nc', inum )
            CALL iom_get  ( inum, jpdom_global, 'aeiu', aeiu, cd_type = 'U', psgn = 1._wp )
            CALL iom_get  ( inum, jpdom_global, 'aeiv', aeiv, cd_type = 'V', psgn = 1._wp )
            CALL iom_close( inum )
            !
         CASE(  30  )                        !--  fixed 3D shape  --!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy induced velocity coef. = F( latitude, longitude, depth )'
            CALL ldf_c2d( 'TRA', zUfac      , inn        , aeiu, aeiv )    ! surface value proportional to scale factor^inn
            CALL ldf_c1d( 'TRA', aeiu(:,:,1), aeiv(:,:,1), aeiu, aeiv )    ! reduction with depth
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_tra_init: wrong choice for nn_aei_ijk_t, the type of space-time variation of aei')
         END SELECT
         !
         IF( .NOT.l_ldfeiv_time ) THEN             !* mask if No time variation
            DO jk = 1, jpkm1
               aeiu(:,:,jk) = aeiu(:,:,jk) * umask(:,:,jk)
               aeiv(:,:,jk) = aeiv(:,:,jk) * vmask(:,:,jk)
            END DO
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_eiv_init


   SUBROUTINE ldf_eiv( kt, paei0, paeiu, paeiv, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!              growth rate of baroclinic instability.
      !!
      !! ** Method  :   coefficient function of the growth rate of baroclinic instability
      !!
      !! Reference : Treguier et al. JPO 1997   ; Held and Larichev JAS 1996
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt             ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm            ! ocean time level indices
      REAL(wp)                        , INTENT(in   ) ::   paei0          ! max value            [m2/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   paeiu, paeiv   ! eiv coefficient      [m2/s]
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zfw, ze3w, zn2, z1_f20, zzaei    ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zn, zah, zhw, zRo, zaeiw   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      zn (:,:) = 0._wp        ! Local initialization
      zhw(:,:) = 5._wp
      zah(:,:) = 0._wp
      zRo(:,:) = 0._wp
      !                       ! Compute lateral diffusive coefficient at T-point
      IF( ln_traldf_triad ) THEN
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
            ! Take the max of N^2 and zero then take the vertical sum
            ! of the square root of the resulting N^2 ( required to compute
            ! internal Rossby radius Ro = .5 * sum_jpk(N) / f
            zn2 = MAX( rn2b(ji,jj,jk), 0._wp )
            zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * e3w(ji,jj,jk,Kmm)
            ! Compute elements required for the inverse time scale of baroclinic
            ! eddies using the isopycnal slopes calculated in ldfslp.F :
            ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
            ze3w = e3w(ji,jj,jk,Kmm) * wmask(ji,jj,jk)
            zah(ji,jj) = zah(ji,jj) + zn2 * wslp2(ji,jj,jk) * ze3w
            zhw(ji,jj) = zhw(ji,jj) + ze3w
         END_3D
      ELSE
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
            ! Take the max of N^2 and zero then take the vertical sum
            ! of the square root of the resulting N^2 ( required to compute
            ! internal Rossby radius Ro = .5 * sum_jpk(N) / f
            zn2 = MAX( rn2b(ji,jj,jk), 0._wp )
            zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * e3w(ji,jj,jk,Kmm)
            ! Compute elements required for the inverse time scale of baroclinic
            ! eddies using the isopycnal slopes calculated in ldfslp.F :
            ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
            ze3w = e3w(ji,jj,jk,Kmm) * wmask(ji,jj,jk)
            zah(ji,jj) = zah(ji,jj) + zn2 * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)   &
               &                            + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) ) * ze3w
            zhw(ji,jj) = zhw(ji,jj) + ze3w
         END_3D
      ENDIF

      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
         zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
         ! Rossby radius at w-point taken betwenn 2 km and  40km
         zRo(ji,jj) = MAX(  2.e3 , MIN( .4 * zn(ji,jj) / zfw, 40.e3 )  )
         ! Compute aeiw by multiplying Ro^2 and T^-1
         zaeiw(ji,jj) = zRo(ji,jj) * zRo(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * tmask(ji,jj,1)
      END_2D

      !                                         !==  Bound on eiv coeff.  ==!
      z1_f20 = 1._wp / (  2._wp * omega * sin( rad * 20._wp )  )
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
         zzaei = MIN( 1._wp, ABS( ff_t(ji,jj) * z1_f20 ) ) * zaeiw(ji,jj)     ! tropical decrease
         zaeiw(ji,jj) = MIN( zzaei , paei0 )                                  ! Max value = paei0
      END_2D
      IF( nn_hls == 1 )   CALL lbc_lnk( 'ldftra', zaeiw(:,:), 'W', 1.0_wp )   ! lateral boundary condition
      !
      DO_2D( 0, 0, 0, 0 )
         paeiu(ji,jj,1) = 0.5_wp * ( zaeiw(ji,jj) + zaeiw(ji+1,jj  ) ) * umask(ji,jj,1)
         paeiv(ji,jj,1) = 0.5_wp * ( zaeiw(ji,jj) + zaeiw(ji  ,jj+1) ) * vmask(ji,jj,1)
      END_2D
      CALL lbc_lnk( 'ldftra', paeiu(:,:,1), 'U', 1.0_wp , paeiv(:,:,1), 'V', 1.0_wp )      ! lateral boundary condition

      DO jk = 2, jpkm1                          !==  deeper values equal the surface one  ==!
         paeiu(:,:,jk) = paeiu(:,:,1) * umask(:,:,jk)
         paeiv(:,:,jk) = paeiv(:,:,1) * vmask(:,:,jk)
      END DO
      !
   END SUBROUTINE ldf_eiv


   SUBROUTINE ldf_eiv_trp( kt, kit000, pu, pv, pw, cdtype, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv_trp  ***
      !!
      !! ** Purpose :   add to the input ocean transport the contribution of
      !!              the eddy induced velocity parametrization.
      !!
      !! ** Method  :   The eddy induced transport is computed from a flux stream-
      !!              function which depends on the slope of iso-neutral surfaces
      !!              (see ldf_slp). For example, in the i-k plan :
      !!                   psi_uw = mk(aeiu) e2u mi(wslpi)   [in m3/s]
      !!                   Utr_eiv = - dk[psi_uw]
      !!                   Vtr_eiv = + di[psi_uw]
      !!                ln_ldfeiv_dia = T : output the associated streamfunction,
      !!                                    velocity and heat transport (call ldf_eiv_dia)
      !!
      !! ** Action  : pu, pv increased by the eiv transport
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kt        ! ocean time-step index
      INTEGER                     , INTENT(in   ) ::   kit000    ! first time step index
      INTEGER                     , INTENT(in   ) ::   Kmm, Krhs ! ocean time level indices
      CHARACTER(len=3)            , INTENT(in   ) ::   cdtype    ! =TRA or TRC (tracer indicator)
      ! TEMP: [tiling] Can be A2D(nn_hls) after all lbc_lnks removed in the nn_hls = 2 case in tra_adv_fct
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu        ! in : 3 ocean transport components   [m3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pv        ! out: 3 ocean transport components   [m3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pw        ! increased by the eiv                [m3/s]
      !!
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zuwk, zuwk1, zuwi, zuwi1   ! local scalars
      REAL(wp) ::   zvwk, zvwk1, zvwj, zvwj1   !   -      -
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zpsi_uw, zpsi_vw
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == kit000 )  THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'ldf_eiv_trp : eddy induced advection on ', cdtype,' :'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   add to velocity fields the eiv component'
         ENDIF
      ENDIF


      zpsi_uw(:,:, 1 ) = 0._wp   ;   zpsi_vw(:,:, 1 ) = 0._wp
      zpsi_uw(:,:,jpk) = 0._wp   ;   zpsi_vw(:,:,jpk) = 0._wp
      !
      DO_3D( nn_hls, nn_hls-1, nn_hls, nn_hls-1, 2, jpkm1 )
         zpsi_uw(ji,jj,jk) = - r1_4 * e2u(ji,jj) * ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk) )   &
            &                                    * ( aeiu (ji,jj,jk-1) + aeiu (ji  ,jj,jk) ) * wumask(ji,jj,jk)
         zpsi_vw(ji,jj,jk) = - r1_4 * e1v(ji,jj) * ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk) )   &
            &                                    * ( aeiv (ji,jj,jk-1) + aeiv (ji,jj  ,jk) ) * wvmask(ji,jj,jk)
      END_3D
      !
      DO_3D_OVR( nn_hls, nn_hls-1, nn_hls, nn_hls-1, 1, jpkm1 )
         pu(ji,jj,jk) = pu(ji,jj,jk) - ( zpsi_uw(ji,jj,jk) - zpsi_uw(ji,jj,jk+1) )
         pv(ji,jj,jk) = pv(ji,jj,jk) - ( zpsi_vw(ji,jj,jk) - zpsi_vw(ji,jj,jk+1) )
      END_3D
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
         pw(ji,jj,jk) = pw(ji,jj,jk) + (  zpsi_uw(ji,jj,jk) - zpsi_uw(ji-1,jj  ,jk)   &
            &                           + zpsi_vw(ji,jj,jk) - zpsi_vw(ji  ,jj-1,jk) )
      END_3D
      !
      !                              ! diagnose the eddy induced velocity and associated heat transport
      IF( ln_ldfeiv_dia .AND. cdtype == 'TRA' )   CALL ldf_eiv_dia( zpsi_uw, zpsi_vw, Kmm )
      !
    END SUBROUTINE ldf_eiv_trp


   SUBROUTINE ldf_eiv_dia( psi_uw, psi_vw, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv_dia  ***
      !!
      !! ** Purpose :   diagnose the eddy induced velocity and its associated
      !!              vertically integrated heat transport.
      !!
      !! ** Method :
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in) ::   psi_uw, psi_vw   ! streamfunction   [m3/s]
      INTEGER                             , INTENT(in) ::   Kmm              ! ocean time level indices
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zztmp   ! local scalar
      REAL(wp), DIMENSION(A2D(nn_hls))     ::   zw2d   ! 2D workspace
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zw3d   ! 3D workspace
      !!----------------------------------------------------------------------
      !
!!gm I don't like this routine....   Crazy  way of doing things, not optimal at all...
!!gm     to be redesigned....
      !                                                  !==  eiv stream function: output  ==!
!!gm      CALL iom_put( "psi_eiv_uw", psi_uw )                 ! output
!!gm      CALL iom_put( "psi_eiv_vw", psi_vw )
      !
      !                                                  !==  eiv velocities: calculate and output  ==!
      !
      zw3d(:,:,jpk) = 0._wp                                    ! bottom value always 0
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )                                  ! e2u e3u u_eiv = -dk[psi_uw]
         zw3d(ji,jj,jk) = ( psi_uw(ji,jj,jk+1) - psi_uw(ji,jj,jk) ) / ( e2u(ji,jj) * e3u(ji,jj,jk,Kmm) )
      END_3D
      CALL iom_put( "uoce_eiv", zw3d )
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )                                  ! e1v e3v v_eiv = -dk[psi_vw]
         zw3d(ji,jj,jk) = ( psi_vw(ji,jj,jk+1) - psi_vw(ji,jj,jk) ) / ( e1v(ji,jj) * e3v(ji,jj,jk,Kmm) )
      END_3D
      CALL iom_put( "voce_eiv", zw3d )
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )                            ! e1 e2 w_eiv = dk[psix] + dk[psix]
         zw3d(ji,jj,jk) = (  psi_vw(ji,jj,jk) - psi_vw(ji  ,jj-1,jk)    &
            &              + psi_uw(ji,jj,jk) - psi_uw(ji-1,jj  ,jk)  ) / e1e2t(ji,jj)
      END_3D
      CALL iom_put( "woce_eiv", zw3d )
      !
      IF( iom_use('weiv_masstr') ) THEN   ! vertical mass transport & its square value
         DO_2D( 0, 0, 0, 0 )
            zw2d(ji,jj) = rho0 * e1e2t(ji,jj)
         END_2D
         DO jk = 1, jpk
            zw3d(:,:,jk) = zw3d(:,:,jk) * zw2d(:,:)
         END DO
         CALL iom_put( "weiv_masstr" , zw3d )
      ENDIF
      !
      IF( iom_use('ueiv_masstr') ) THEN
         zw3d(:,:,:) = 0.e0
         DO jk = 1, jpkm1
            zw3d(:,:,jk) = rho0 * ( psi_uw(:,:,jk+1) - psi_uw(:,:,jk) )
         END DO
         CALL iom_put( "ueiv_masstr", zw3d )                  ! mass transport in i-direction
      ENDIF
      !
      zztmp = 0.5_wp * rho0 * rcp
      IF( iom_use('ueiv_heattr') .OR. iom_use('ueiv_heattr3d') ) THEN
        zw2d(:,:)   = 0._wp
        zw3d(:,:,:) = 0._wp
        DO_3D( 0, 0, 0, 0, 1, jpkm1 )
           zw3d(ji,jj,jk) = zw3d(ji,jj,jk) + ( psi_uw(ji,jj,jk+1)          - psi_uw(ji  ,jj,jk)            )   &
              &                            * ( ts    (ji,jj,jk,jp_tem,Kmm) + ts    (ji+1,jj,jk,jp_tem,Kmm) )
           zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
        END_3D
        CALL iom_put( "ueiv_heattr"  , zztmp * zw2d )                  ! heat transport in i-direction
        CALL iom_put( "ueiv_heattr3d", zztmp * zw3d )                  ! heat transport in i-direction
      ENDIF
      !
      IF( iom_use('veiv_masstr') ) THEN
         zw3d(:,:,:) = 0.e0
         DO jk = 1, jpkm1
            zw3d(:,:,jk) = rho0 * ( psi_vw(:,:,jk+1) - psi_vw(:,:,jk) )
         END DO
         CALL iom_put( "veiv_masstr", zw3d )                  ! mass transport in i-direction
      ENDIF
      !
      zw2d(:,:)   = 0._wp
      zw3d(:,:,:) = 0._wp
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zw3d(ji,jj,jk) = zw3d(ji,jj,jk) + ( psi_vw(ji,jj,jk+1)          - psi_vw(ji,jj  ,jk)            )   &
            &                            * ( ts    (ji,jj,jk,jp_tem,Kmm) + ts    (ji,jj+1,jk,jp_tem,Kmm) )
         zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
      END_3D
      CALL iom_put( "veiv_heattr"  , zztmp * zw2d )                  !  heat transport in j-direction
      CALL iom_put( "veiv_heattr3d", zztmp * zw3d )                  !  heat transport in j-direction
      !
      IF( iom_use( 'sophteiv' ) )   CALL dia_ptr_hst( jp_tem, 'eiv', 0.5 * zw3d )
      !
      zztmp = 0.5_wp * 0.5
      IF( iom_use('ueiv_salttr') .OR. iom_use('ueiv_salttr3d')) THEN
        zw2d(:,:) = 0._wp
        zw3d(:,:,:) = 0._wp
        DO_3D( 0, 0, 0, 0, 1, jpkm1 )
           zw3d(ji,jj,jk) = zw3d(ji,jj,jk) * ( psi_uw(ji,jj,jk+1)          - psi_uw(ji  ,jj,jk)            )   &
              &                            * ( ts    (ji,jj,jk,jp_sal,Kmm) + ts    (ji+1,jj,jk,jp_sal,Kmm) )
           zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
        END_3D
        CALL iom_put( "ueiv_salttr", zztmp * zw2d )                  ! salt transport in i-direction
        CALL iom_put( "ueiv_salttr3d", zztmp * zw3d )                ! salt transport in i-direction
      ENDIF
      zw2d(:,:) = 0._wp
      zw3d(:,:,:) = 0._wp
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zw3d(ji,jj,jk) = zw3d(ji,jj,jk) + ( psi_vw(ji,jj,jk+1)          - psi_vw(ji,jj  ,jk)            )   &
            &                            * ( ts    (ji,jj,jk,jp_sal,Kmm) + ts    (ji,jj+1,jk,jp_sal,Kmm) )
         zw2d(ji,jj) = zw2d(ji,jj) + zw3d(ji,jj,jk)
      END_3D
      CALL iom_put( "veiv_salttr"  , zztmp * zw2d )                  !  salt transport in j-direction
      CALL iom_put( "veiv_salttr3d", zztmp * zw3d )                  !  salt transport in j-direction
      !
      IF( iom_use( 'sopsteiv' ) ) CALL dia_ptr_hst( jp_sal, 'eiv', 0.5 * zw3d )
      !
      !
   END SUBROUTINE ldf_eiv_dia

   !!======================================================================
END MODULE ldftra
