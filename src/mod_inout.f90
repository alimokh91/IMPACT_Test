!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* Apr 2012                                                                                                  *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch                    *
!* October 2014                                                                                              *
!*************************************************************************************************************

!> module containing subroutines for reading and writing data to hdf5 files. It uses modules mod_dims, mod_vars,
!! mod_exchange, mod_lib and HDF5.
!! @todo: information on the immersed boundary could be written to files with routines that are added here.
!!        If not, there should be routines in the structural solver.
MODULE mod_inout
 
 
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib
  USE HDF5
  USE MPI
 
  PRIVATE
 
  PUBLIC write_fields
  PUBLIC write_restart, read_restart
  PUBLIC write_hdf, read_hdf, read2_hdf, write_hdf_velall
  PUBLIC write_2D_hdf, read2_2D_hdf
  PUBLIC write_1D_hdf, read2_1D_hdf
  PUBLIC write_stats_hdf_4D, read_stats_hdf_4D
  PUBLIC write_stats_hdf_2D, read_stats_hdf_2D
 
  PUBLIC write_hdf_infoREAL, read_hdf_infoREAL
  PUBLIC write_hdf_infoINT , read_hdf_infoINT
  PUBLIC write_hdf_infoLOG , read_hdf_infoLOG
  PUBLIC start_hdf5_for_testing, stop_hdf5_for_testing !bbecsek

  PUBLIC filespace_props
  PUBLIC lambda
  PUBLIC vorticity_2D
 
  PUBLIC write_xdmf_xml
  PUBLIC write_xdmf_timecollection
  PUBLIC write_2D_xdmf, write_3D_xdmf
 
  INTEGER(HID_T)                 ::  file_id, plist_id, dset_id
  INTEGER(HID_T)                 ::  filespace, memspace
  INTEGER(HID_T)                 ::  memtypeREAL, memtypeINT
 
  INTEGER(HSIZE_T )              ::  dims_file  (1:3) ! wird global eingefuehrt, um Probleme bei Uebergabe als Argument nach "CALL filespace_props" zu vermeiden.
  INTEGER(HSSIZE_T)              ::  offset_file(1:3)
 
  CONTAINS
 
!!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
 
 

  SUBROUTINE write_fields
 
  IMPLICIT NONE
 
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  REAL :: write_field(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1)
  
  CHARACTER(LEN=8)       ::  count_char
  CHARACTER(LEN=2)       ::  id,phs
  CHARACTER(LEN=1)       ::  conc_number
  CHARACTER(LEN=50)      ::  write_dir

  IF (rank == 0) WRITE(*,'(a)') 'writing fields ...'
  
  !===========================================================================================================
  !=== Ausschrieb-Nr. als String fuer File-Namen =============================================================
  !===========================================================================================================
  CALL num_to_string(8,write_count,count_char)
  !===========================================================================================================
  
  
  CALL exchange(1,1,vel(b1L,b2L,b3L,1))
  CALL exchange(2,2,vel(b1L,b2L,b3L,2))
  CALL exchange(3,3,vel(b1L,b2L,b3L,3))
  IF (write_force_yes) THEN
    CALL exchange(1,1,fd(b1L,b2L,b3L,1))
    CALL exchange(2,2,fd(b1L,b2L,b3L,2))
    CALL exchange(3,3,fd(b1L,b2L,b3L,3))
  END IF
  
  !===========================================================================================================
  !=== Interpolieren / Schreiben =============================================================================
  !===========================================================================================================
  IF (dimens == 3) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1) ! TEST!!! Verallgemeinerte Interpolation verwenden? Mit compute_stats teilen?
              END DO
           END DO
        END DO
     END DO
     IF (write_large) CALL write_hdf('velX_large_'//count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,1)) ! TEST!!! velX --> vel1, etc. umbenennen!!
     IF (write_med  ) CALL write_hdf('velX_med_'  //count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,1))
     IF (write_small) CALL write_hdf('velX_small_'//count_char,'velX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2L,b3L,1))
  ELSE
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              bc33(i,j,1) = bc33(i,j,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
     CALL write_2D_hdf('velX_'//count_char,'velX',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
              END DO
           END DO
        END DO
     END DO
     IF (write_large) CALL write_hdf('velY_large_'//count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,2))
     IF (write_med  ) CALL write_hdf('velY_med_'  //count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,2))
     IF (write_small) CALL write_hdf('velY_small_'//count_char,'velY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2L,b3L,2))
  ELSE
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              bc33(i,j,2) = bc33(i,j,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
     CALL write_2D_hdf('velY_'//count_char,'velY',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,2))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,3) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 nl(i,j,k,3) = nl(i,j,k,3) + cIwp(kk,k)*vel(i,j,k+kk,3)
              END DO
           END DO
        END DO
     END DO
     IF (write_large) CALL write_hdf('velZ_large_'//count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,3))
     IF (write_med  ) CALL write_hdf('velZ_med_'  //count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,3))
     IF (write_small) CALL write_hdf('velZ_small_'//count_char,'velZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2L,b3L,3))
  END IF
  !===========================================================================================================
  IF (1 == 2 .AND. dimens == 3) CALL write_hdf_velall('velA_'//count_char,S1p,S2p,S3p,N1p,N2p,N3p,0,nl)
  !===========================================================================================================
  IF (write_lambda2_yes .AND. dimens == 3) THEN
     CALL lambda
     IF (write_large) CALL write_hdf('lamb_large_'//count_char,'lamb',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,res)
     IF (write_med  ) CALL write_hdf('lamb_med_'  //count_char,'lamb',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,res)
     IF (write_small) CALL write_hdf('lamb_small_'//count_char,'lamb',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,res)
  ELSE IF (write_lambda2_yes) THEN
     CALL vorticity_2D
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,1) = res(i,j,k)
        END DO
     END DO
     CALL write_2D_hdf('vort_'//count_char,'vort',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
  END IF
  !===========================================================================================================
  IF (dimens == 3) THEN
     IF (write_large) CALL write_hdf('pre_large_'//count_char,'pre',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,pre)
     IF (write_med  ) CALL write_hdf('pre_med_'  //count_char,'pre',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,pre)
     IF (write_small) CALL write_hdf('pre_small_'//count_char,'pre',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,pre)
  ELSE
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,1) = pre(i,j,k)
        END DO
     END DO
     CALL write_2D_hdf('pre_'//count_char,'pre',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
  END IF
  !===========================================================================================================
  IF (dimens == 3 .AND. write_force_yes) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,1) = cIup(d1L,i)*fd(i+d1L,j,k,1)
!pgi$ unroll = n:8
              DO ii = d1L+1, d1U
                 nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*fd(i+ii,j,k,1) ! TEST!!! Verallgemeinerte Interpolation verwenden? Mit compute_stats teilen?
              END DO
           END DO
        END DO
     END DO
    IF (write_large) CALL write_hdf('forceX_large_'//count_char,'forceX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,1))
    IF (write_med  ) CALL write_hdf('forceX_med_'  //count_char,'forceX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,1))
    IF (write_small) CALL write_hdf('forceX_small_'//count_char,'forceX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2l,b3L,1))
  ELSE IF (write_force_yes) THEN
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,1) = cIup(d1L,i)*fd(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              bc33(i,j,1) = bc33(i,j,1) + cIup(ii,i)*fd(i+ii,j,k,1)
           END DO
        END DO
     END DO
    !k = 1
    !DO j = S2p, N2p
    !  DO i = S1p, N1p
    !    bc33(i,j,1) = fd(i,j,k,1)
    !  END DO
    !END DO
    CALL write_2D_hdf('forceX_'//count_char,'forceX',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,1))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3 .AND. write_force_yes) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,2) = cIvp(d2L,j)*fd(i,j+d2L,k,2)
!pgi$ unroll = n:8
              DO jj = d2L+1, d2U
                 nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*fd(i,j+jj,k,2)
              END DO
           END DO
        END DO
     END DO
    IF (write_large) CALL write_hdf('forceY_large_'//count_char,'forceY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,2))
    IF (write_med  ) CALL write_hdf('forceY_med_'  //count_char,'forceY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,2))
    IF (write_small) CALL write_hdf('forceY_small_'//count_char,'forceY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2l,b3L,2))
  ELSE IF (write_force_yes) THEN
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           bc33(i,j,2) = cIvp(d2L,j)*fd(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              bc33(i,j,2) = bc33(i,j,2) + cIvp(jj,j)*fd(i,j+jj,k,2)
           END DO
        END DO
     END DO
    !k = 1
    !DO j = S2p, N2p
    !  DO i = S1p, N1p
    !    bc33(i,j,2) = fd(i,j,k,2)
    !  END DO
    !END DO
    CALL write_2D_hdf('forceY_'//count_char,'forceY',N1,N2,S1p,S2p,N1p,N2p,iShift,jShift,-3,bc33(1,1,2))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3 .AND. write_force_yes) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              nl(i,j,k,3) = cIwp(d3L,k)*fd(i,j,k+d3L,3)
!pgi$ unroll = n:8
              DO kk = d3L+1, d3U
                 nl(i,j,k,3) = nl(i,j,k,3) + cIwp(kk,k)*fd(i,j,k+kk,3)
              END DO
           END DO
        END DO
     END DO
    IF (write_large) CALL write_hdf('forceZ_large_'//count_char,'forceZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,nl(b1L,b2L,b3L,3))
    IF (write_med  ) CALL write_hdf('forceZ_med_'  //count_char,'forceZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_med  ,nl(b1L,b2L,b3L,3))
    IF (write_small) CALL write_hdf('forceZ_small_'//count_char,'forceZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_small,nl(b1L,b2l,b3L,3))
  END IF
  !===========================================================================================================

  !IF (write_covariance_yes ) THEN

  !  if (dtime_out_kalm.ne.0.0) then
  !     phase = mod(write_kalm_count-1,intervals) + 1
  !     CALL num_to_string(2,phase,phs)
  !     write_dir = './kf_result/phase_'//phs//'/'
  !     CALL num_to_string(8,write_kalm_count,count_char)

  !     write_field = 0.
  !     if (associated(kalman_first)) then
  !       DO kk = (lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1, &
  !                ubound(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3)+1
  !         DO jj = (lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1, &
  !                  ubound(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1
  !           DO ii = (lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1, &
  !                    ubound(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1

  !             i =      ii-((lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1)
  !             i = i + (jj-((lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1))* &
  !                     (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1)
  !             i = i + (kk-((lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1))* &
  !                     (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1) * &
  !                     (size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)

  !             write_field(ii,jj,kk,1) = kalman_first%K(3*i+1,1)
  !           END DO
  !         END DO
  !       END DO
  !     end if
  !     CALL write_hdf(trim(write_dir)//'gainX_phase'//phs//'_'//count_char,'gainX',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_field(b1L,b2L,b3L,1))
  !     write_field = 0.
  !     if (associated(kalman_first)) then
  !       DO kk = (lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1, &
  !                ubound(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3)+1
  !         DO jj = (lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1, &
  !                  ubound(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1
  !           DO ii = (lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1, &
  !                    ubound(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1

  !             i =      ii-((lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1)
  !             i = i + (jj-((lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1))* &
  !                     (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1)
  !             i = i + (kk-((lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1))* &
  !                     (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1) * &
  !                     (size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)

  !             write_field(ii,jj,kk,1) = kalman_first%K(3*i+2,1)
  !           END DO
  !         END DO
  !       END DO
  !     end if
  !     CALL write_hdf(trim(write_dir)//'gainY_phase'//phs//'_'//count_char,'gainY',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_field(b1L,b2L,b3L,1))
  !     write_field = 0.
  !     if (associated(kalman_first)) then
  !       DO kk = (lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1, &
  !                ubound(mri_inst%mri%velocity_mean%array,5)*kalman_num_spatial_refinements(3)+1
  !         DO jj = (lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1, &
  !                  ubound(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1
  !           DO ii = (lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1, &
  !                    ubound(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1

  !             i =      ii-((lbound(mri_inst%mri%velocity_mean%array,3)-1)*kalman_num_spatial_refinements(1)+1)
  !             i = i + (jj-((lbound(mri_inst%mri%velocity_mean%array,4)-1)*kalman_num_spatial_refinements(2)+1))* &
  !                     (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1)
  !             i = i + (kk-((lbound(mri_inst%mri%velocity_mean%array,5)-1)*kalman_num_spatial_refinements(3)+1))* &
  !                     (size(mri_inst%mri%velocity_mean%array,3)*kalman_num_spatial_refinements(1)+1) * &
  !                     (size(mri_inst%mri%velocity_mean%array,4)*kalman_num_spatial_refinements(2)+1)

  !             write_field(ii,jj,kk,1) = kalman_first%K(3*i+3,1)
  !           END DO
  !         END DO
  !       END DO
  !     end if
  !     CALL write_hdf(trim(write_dir)//'gainZ_phase'//phs//'_'//count_char,'gainZ',S1p,S2p,S3p,N1p,N2p,N3p,0,stride_large,write_field(b1L,b2L,b3L,1))

  !  end if
 
  !END IF

  write_count    = write_count   + 1
  time_out_vect  = time_out_vect + dtime_out_vect
  write_out_vect = .FALSE.
 
 
  END SUBROUTINE write_fields
 
 
 
 
 
 
 
 
 
 
 
  SUBROUTINE write_restart
  
  IMPLICIT NONE
  
  CHARACTER(LEN=3)       ::  next_restart_char
  
  
  IF (write_restart_yes) THEN
     
     IF (rank == 0) WRITE(*,'(a,i3,a)') 'writing data for restart',restart,' ...'
     
     !========================================================================================================
     !=== neue Restart-Nr. als String fuer File-Namen ========================================================
     !========================================================================================================
     CALL num_to_string(3,restart,next_restart_char)
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Schreiben ==========================================================================================
     !========================================================================================================
                      CALL write_hdf('velX_restart'//next_restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,(/1,1,1/),vel(b1L,b2L,b3L,1))
                      CALL write_hdf('velY_restart'//next_restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,(/1,1,1/),vel(b1L,b2L,b3L,2))
     IF (dimens == 3) CALL write_hdf('velZ_restart'//next_restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,(/1,1,1/),vel(b1L,b2L,b3L,3))
     !--------------------------------------------------------------------------------------------------------
     CALL write_hdf('pre_restart'//next_restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),pre)
     !--------------------------------------------------------------------------------------------------------
     !========================================================================================================
     
  END IF
  
  
  END SUBROUTINE write_restart
  
  
  
  
  
  
  
  
  
  
  !> subroutine for reading restart files.
  SUBROUTINE read_restart
  
  IMPLICIT NONE
  
  REAL                   ::  old_dtime_out_kalm
  REAL                   ::  old_dtime_out_scal
  REAL                   ::  old_dtime_out_vect

  CHARACTER(LEN=1)       ::  conc_number
  
  
  IF (rank == 0) WRITE(*,'(a,i3,a)') 'reading data for restart',restart,' ...'
  IF (rank == 0) WRITE(*,*)
  
 
                   CALL read2_hdf('velX_restart'//restart_char,'velX_restart',S11B,S21B,S31B,N11B,N21B,N31B,1,vel(b1L,b2L,b3L,1))
                   CALL read2_hdf('velY_restart'//restart_char,'velY_restart',S12B,S22B,S32B,N12B,N22B,N32B,2,vel(b1L,b2L,b3L,2))
  IF (dimens == 3) CALL read2_hdf('velZ_restart'//restart_char,'velZ_restart',S13B,S23B,S33B,N13B,N23B,N33B,3,vel(b1L,b2L,b3L,3))
  !-----------------------------------------------------------------------------------------------------------
  CALL read2_hdf('pre_restart'//restart_char,'pre_restart',S1p,S2p,S3p,N1p,N2p,N3p,0,pre)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------

         ! --- mjohn 070213 - allow for possibility of adaptation of dtime_out_vect or dtime_out_scal by config.txt after restart
         ! compute next time_out_vect and next time_out_scal by subtracting last dtime (stored in restart file) and adding new dtime (read from config.txt)

         ! Setup file access property list with parallel I/O access
         CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
         CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

         ! Open file collectively
         CALL h5fopen_f('velX_restart'//restart_char//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
         CALL h5pclose_f(plist_id,herror)

         CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
         CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
         !-----------------------------------------------------------------------------------------------------------
         ! Open the dataset:
         CALL h5dopen_f(file_id,'velX_restart',dset_id,herror)
         CALL read_hdf_infoREAL(1,.TRUE. ,.TRUE.,'dtime_out_kalm' ,scalar=old_dtime_out_kalm)
         CALL read_hdf_infoREAL(1,.TRUE. ,.TRUE.,'dtime_out_scal' ,scalar=old_dtime_out_scal)
         CALL read_hdf_infoREAL(1,.TRUE. ,.TRUE.,'dtime_out_vect' ,scalar=old_dtime_out_vect)
         CALL h5dclose_f(dset_id,herror)
         ! Close the file
         CALL h5fclose_f(file_id,herror) 
         !===========================================================================================================

         time_out_vect = (time_out_vect-old_dtime_out_vect) + dtime_out_vect
         time_out_scal = (time_out_scal-old_dtime_out_scal) + dtime_out_scal
         time_out_kalm = (time_out_kalm-old_dtime_out_kalm) + dtime_out_kalm

         ! catch values less than or equal to 'time'
         if(time_out_vect .le. time) then
            time_out_vect = time + dtime
         end if
         if(time_out_scal .le. time) then
            time_out_scal = time + dtime
         end if
         if(time_out_kalm .le. time) then
            time_out_kalm = time + dtime
         end if

  
  END SUBROUTINE read_restart
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,stride,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3
  
  INTEGER, INTENT(in)    ::  NN1
  INTEGER, INTENT(in)    ::  NN2
  INTEGER, INTENT(in)    ::  NN3
  
  INTEGER, INTENT(in)    ::  vel_dir
  INTEGER, INTENT(in)    ::  stride(1:3)
  
  REAL   , INTENT(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER(HSIZE_T )      ::  stride_mem(1:3)
  INTEGER(HSIZE_T )      ::  dims_mem  (1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)      ::  offset_mem(1:3)
  
  INTEGER                ::  S1w, M1w, dim1
  INTEGER                ::  S2w, M2w, dim2
  INTEGER                ::  S3w, M3w, dim3
  
  LOGICAL                ::  attr_yes
  
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! Nur fuer Druckgitter! (sonst wird es richtig kompliziert, siehe Multigrid fuer Helmholtz-Problem!).
  IF (vel_dir /= 0) THEN
     IF (stride(vel_dir) /= 1) THEN
        IF (rank == 0) WRITE(*,*) 'ERROR! Cannot write velocity field with stride /= 1 in the corresponding direction!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
  END IF
  IF ((MOD((N1-1),stride(1)) /= 0) .OR. (MOD((N2-1),stride(2)) /= 0) .OR. (MOD((N3-1),stride(3)) /= 0)) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot write field with this stride!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  CALL filespace_props(vel_dir,stride,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  dims_data  = (/((NN1-SS1)/stride(1)+1),((NN2-SS2)/stride(2)+1),((NN3-SS3)/stride(3)+1)/)
  offset_mem = (/(SS1-b1L),(SS2-b2L),(SS3-b3L)/)
  stride_mem = (/stride(1),stride(2),stride(3)/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_kalm'   ,scalar=dtime_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  !--- ddemarinis: DA additions
  CALL write_hdf_infoINT (1         ,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL write_hdf_infoREAL(1         ,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL write_hdf_infoREAL(1         ,.TRUE. ,attr_yes,'dtime_out_kalm'   ,scalar=dtime_out_kalm   )
  CALL write_hdf_infoLOG (1         ,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1  M2  M3 '      ,array =(/M1 ,M2 ,M3 /)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB1,NB2,NB3/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls1,ls2,ls3/)  )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'CFL   '           ,scalar=CFL              )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'thetaL'           ,scalar=thetaL           )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'epsU  '           ,scalar=epsU             )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
  
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iL',array=(/BC_1L_global,BC_2L_global,BC_3L_global/))
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iU',array=(/BC_1U_global,BC_2U_global,BC_3U_global/))
  
  CALL write_hdf_infoINT (dim_ncb1c,.FALSE.,attr_yes,'ncb1c',array=ncb1c    )
  CALL write_hdf_infoINT (dim_ncb1g,.FALSE.,attr_yes,'ncb1g',array=ncb1g    )
  CALL write_hdf_infoINT (dim_ncb1d,.FALSE.,attr_yes,'ncb1d',array=ncb1d    )
  
  CALL write_hdf_infoINT (dim_ncb2c,.FALSE.,attr_yes,'ncb2c',array=ncb2c    )
  CALL write_hdf_infoINT (dim_ncb2g,.FALSE.,attr_yes,'ncb2g',array=ncb2g    )
  CALL write_hdf_infoINT (dim_ncb2d,.FALSE.,attr_yes,'ncb2d',array=ncb2d    )
  
  CALL write_hdf_infoINT (dim_ncb3c,.FALSE.,attr_yes,'ncb3c',array=ncb3c    )
  CALL write_hdf_infoINT (dim_ncb3g,.FALSE.,attr_yes,'ncb3g',array=ncb3g    )
  CALL write_hdf_infoINT (dim_ncb3d,.FALSE.,attr_yes,'ncb3d',array=ncb3d    )
  
  CALL write_hdf_infoREAL( M1      ,.FALSE.,attr_yes,'y1p'  ,array=y1p(1:M1))
  CALL write_hdf_infoREAL( M2      ,.FALSE.,attr_yes,'y2p'  ,array=y2p(1:M2))
  CALL write_hdf_infoREAL( M3      ,.FALSE.,attr_yes,'y3p'  ,array=y3p(1:M3))
  
  CALL write_hdf_infoREAL((M1+1)   ,.FALSE.,attr_yes,'y1u'  ,array=y1u(0:M1))
  CALL write_hdf_infoREAL((M2+1)   ,.FALSE.,attr_yes,'y2v'  ,array=y2v(0:M2))
  CALL write_hdf_infoREAL((M3+1)   ,.FALSE.,attr_yes,'y3w'  ,array=y3w(0:M3))
  !--- bbecsek: FSI additions
  CALL write_hdf_infoREAL(1        ,.TRUE. ,attr_yes,'U_ref',scalar=U_ref   )
  CALL write_hdf_infoREAL(1        ,.TRUE. ,attr_yes,'L_ref',scalar=L_ref   )
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror,stride_mem)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  IF (vel_dir == 1) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1u(S1w::stride(1)))
  ELSE
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w::stride(1)))
  END IF
  IF (vel_dir == 2) THEN
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2v(S2w::stride(2)))
  ELSE
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w::stride(2)))
  END IF
  IF (vel_dir == 3) THEN
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3w(S3w::stride(3)))
  ELSE
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3p(S3w::stride(3)))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_hdf
  
  
  
  SUBROUTINE write_hdf_velall(filename,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  
  INTEGER, INTENT(in)      ::  SS1
  INTEGER, INTENT(in)      ::  SS2
  INTEGER, INTENT(in)      ::  SS3
  
  INTEGER, INTENT(in)      ::  NN1
  INTEGER, INTENT(in)      ::  NN2
  INTEGER, INTENT(in)      ::  NN3
  
  INTEGER, INTENT(in)      ::  vel_dir
  
  REAL   , INTENT(in)      ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  INTEGER(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)        ::  offset_mem(1:3)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  INTEGER                  ::  S3w, M3w, dim3
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  CALL filespace_props(vel_dir,(/1,1,1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  dims_data  = (/(NN1-SS1+1   ),(NN2-SS2+1   ),(NN3-SS3+1   )/)
  offset_mem = (/(SS1-b1L     ),(SS2-b2L     ),(SS3-b3L     )/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'velX',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_kalm'   ,scalar=dtime_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1  M2  M3 '      ,array =(/M1 ,M2 ,M3 /)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'S1w S2w S3w'      ,array =(/S1w,S2w,S3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'M1w M2w M3w'      ,array =(/M1w,M2w,M3w/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'NB1 NB2 NB3'      ,array =(/NB1,NB2,NB3/)  )
  CALL write_hdf_infoINT (3     ,.FALSE.,attr_yes,'ls1 ls2 ls3'      ,array =(/ls1,ls2,ls3/)  )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'CFL   '           ,scalar=CFL              )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'thetaL'           ,scalar=thetaL           )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'epsU  '           ,scalar=epsU             )
  
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'Euler_yes'        ,scalar=Euler_yes        )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'twostep_yes'      ,scalar=twostep_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'mapping_yes'      ,scalar=mapping_yes      )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'upwind_yes'       ,scalar=upwind_yes       )
  
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iL',array=(/BC_1L_global,BC_2L_global,BC_3L_global/))
  CALL write_hdf_infoINT (3,.FALSE.,attr_yes,'BC_iU',array=(/BC_1U_global,BC_2U_global,BC_3U_global/))
  
  CALL write_hdf_infoINT (dim_ncb1c,.FALSE.,attr_yes,'ncb1c',array=ncb1c    )
  CALL write_hdf_infoINT (dim_ncb1g,.FALSE.,attr_yes,'ncb1g',array=ncb1g    )
  CALL write_hdf_infoINT (dim_ncb1d,.FALSE.,attr_yes,'ncb1d',array=ncb1d    )
  
  CALL write_hdf_infoINT (dim_ncb2c,.FALSE.,attr_yes,'ncb2c',array=ncb2c    )
  CALL write_hdf_infoINT (dim_ncb2g,.FALSE.,attr_yes,'ncb2g',array=ncb2g    )
  CALL write_hdf_infoINT (dim_ncb2d,.FALSE.,attr_yes,'ncb2d',array=ncb2d    )
  
  CALL write_hdf_infoINT (dim_ncb3c,.FALSE.,attr_yes,'ncb3c',array=ncb3c    )
  CALL write_hdf_infoINT (dim_ncb3g,.FALSE.,attr_yes,'ncb3g',array=ncb3g    )
  CALL write_hdf_infoINT (dim_ncb3d,.FALSE.,attr_yes,'ncb3d',array=ncb3d    )
  
  CALL write_hdf_infoREAL( M1      ,.FALSE.,attr_yes,'y1p'  ,array=y1p(1:M1))
  CALL write_hdf_infoREAL( M2      ,.FALSE.,attr_yes,'y2p'  ,array=y2p(1:M2))
  CALL write_hdf_infoREAL( M3      ,.FALSE.,attr_yes,'y3p'  ,array=y3p(1:M3))
  
  CALL write_hdf_infoREAL((M1+1)   ,.FALSE.,attr_yes,'y1u'  ,array=y1u(0:M1))
  CALL write_hdf_infoREAL((M2+1)   ,.FALSE.,attr_yes,'y2v'  ,array=y2v(0:M2))
  CALL write_hdf_infoREAL((M3+1)   ,.FALSE.,attr_yes,'y3w'  ,array=y3w(0:M3))
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,1),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'velY',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,2),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'velZ',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi(b1L,b2L,b3L,3),dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (1 == 1) THEN ! TEST!!! wird u.a. gebraucht, um 'velabs' zu speichern.
  ! Create file space / memory space:
  CALL h5screate_simple_f(3,dims_file,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem ,memspace ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create the dataset:
  CALL h5dcreate_f(file_id,'pre',H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,pre,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (vel_dir == 1) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1u(S1w:M1w))
  ELSE
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
  END IF
  IF (vel_dir == 2) THEN
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2v(S2w:M2w))
  ELSE
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w:M2w))
  END IF
  IF (vel_dir == 3) THEN
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3w(S3w:M3w))
  ELSE
     CALL write_hdf_infoREAL(dim3,.FALSE.,.FALSE.,'VectorZ',array=y3p(S3w:M3w))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_hdf_velall
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! evtl. schoener: write_hdf_2D??
  ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
  SUBROUTINE write_2D_hdf(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  N1
  INTEGER, INTENT(in)      ::  N2
  
  INTEGER, INTENT(in)      ::  SS1
  INTEGER, INTENT(in)      ::  SS2
  
  INTEGER, INTENT(in)      ::  NN1
  INTEGER, INTENT(in)      ::  NN2
  
  INTEGER, INTENT(in)      ::  iShift
  INTEGER, INTENT(in)      ::  jShift
  
  INTEGER, INTENT(in)      ::  dir
  
  REAL   , INTENT(in)      ::  phi(1:N1,1:N2)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
  INTEGER(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  IF (ABS(dir) == 1) THEN
     S1w = 2  + ls2
     S2w = 2  + ls3
     
     M1w = M2 + ls2
     M2w = M3 + ls3
     
     IF (BC_2L_global /= -1) THEN
        S1w = 1
        M1w = M2
     END IF
     IF (BC_3L_global /= -1) THEN
        S2w = 1
        M2w = M3
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (ABS(dir) == 2) THEN
     S1w = 2  + ls1
     S2w = 2  + ls3
     
     M1w = M1 + ls1
     M2w = M3 + ls3
     
     IF (BC_1L_global /= -1) THEN
        S1w = 1
        M1w = M1
     END IF
     IF (BC_3L_global /= -1) THEN
        S2w = 1
        M2w = M3
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (ABS(dir) == 3) THEN
     S1w = 2  + ls1
     S2w = 2  + ls2
     
     M1w = M1 + ls1
     M2w = M2 + ls2
     
     IF (BC_1L_global /= -1) THEN
        S1w = 1
        M1w = M1
     END IF
     IF (BC_2L_global /= -1) THEN
        S2w = 1
        M2w = M2
     END IF
  END IF
  !===========================================================================================================
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  
  dims_file   = (/ dim1 , dim2 /)
  dims_mem    = (/ N1   , N2   /)
  dims_data   = (/(NN1-SS1+1),(NN2-SS2+1)/)
  offset_mem  = (/ SS1-1, SS2-1/)
  offset_file = (/iShift,jShift/)
  
  IF (.NOT. ((dir == -1 .AND. iB(1,1) == 1  ) .OR. (dir == -2 .AND. iB(2,1) == 1  ) .OR. (dir == -3 .AND. iB(3,1) == 1   ) .OR.    &
        &    (dir ==  1 .AND. iB(1,1) == NB1) .OR. (dir ==  2 .AND. iB(2,1) == NB2) .OR. (dir ==  3 .AND. iB(3,1) == NB3))) THEN
     dims_data   = (/0,0/)
     offset_mem  = (/0,0/)
     offset_file = (/0,0/)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(2,dims_file,filespace,herror)
  CALL h5screate_simple_f(2,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_kalm'   ,scalar=dtime_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  
  CALL write_hdf_infoINT (2     ,.FALSE.,attr_yes,'S1w S2w'          ,array =(/S1w,S2w/)      )
  CALL write_hdf_infoINT (2     ,.FALSE.,attr_yes,'M1w M2w'          ,array =(/M1w,M2w/)      )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'Re'               ,scalar=Re               )
  !--- bbecsek: FSI additions
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'U_ref'            ,scalar=U_ref            )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'L_ref'            ,scalar=L_ref            )
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  IF (ABS(dir) == 1) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorY',array=y2p(S1w:M1w))
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorZ',array=y3p(S2w:M2w))
  END IF
  IF (ABS(dir) == 2) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorZ',array=y3p(S2w:M2w))
  END IF
  IF (ABS(dir) == 3) THEN
     CALL write_hdf_infoREAL(dim1,.FALSE.,.FALSE.,'VectorX',array=y1p(S1w:M1w))
     CALL write_hdf_infoREAL(dim2,.FALSE.,.FALSE.,'VectorY',array=y2p(S2w:M2w))
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_2D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_1D_hdf(filename,dsetname,NN,SSp,NNp,iShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  NN
  INTEGER, INTENT(in)      ::  SSp
  INTEGER, INTENT(in)      ::  NNp
  INTEGER, INTENT(in)      ::  iShift
  INTEGER, INTENT(in)      ::  dir
  REAL   , INTENT(in)      ::  phi(1:NN)
  
  INTEGER                  ::  SSw, MMw, dims
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
  INTEGER(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  IF (dir == 1) THEN
     SSw = 2  + ls1
     MMw = M1 + ls1
     IF (BC_1L_global /= -1) THEN
        SSw = 1
        MMw = M1
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dir == 2) THEN
     SSw = 2  + ls2
     MMw = M2 + ls2
     IF (BC_2L_global /= -1) THEN
        SSw = 1
        MMw = M2
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dir == 3) THEN
     SSw = 2  + ls3
     MMw = M3 + ls3
     IF (BC_3L_global /= -1) THEN
        SSw = 1
        MMw = M3
     END IF
  END IF
  !===========================================================================================================
  dims = MMw-SSw+1
  
  dims_file   = (/dims     /)
  dims_mem    = (/NN       /)
  dims_data   = (/NNp-SSp+1/)
  offset_mem  = (/SSp-1    /)
  offset_file = (/iShift   /)
  
  IF (.NOT. ((dir == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) .OR.     &
        &    (dir == 2 .AND. iB(1,1) == 1 .AND. iB(3,1) == 1) .OR.     &
        &    (dir == 3 .AND. iB(1,1) == 1 .AND. iB(2,1) == 1))) THEN
     dims_data   = (/0/)
     offset_mem  = (/0/)
     offset_file = (/0/)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(1,dims_file,filespace,herror)
  CALL h5screate_simple_f(1,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'SSw'              ,scalar=SSw              ) ! TEST!!! SSw, MMw unschoen ...
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'MMw'              ,scalar=MMw              )
  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  IF (dir == 1) CALL write_hdf_infoREAL(dims,.FALSE.,.FALSE.,'VectorX',array=y1p(SSw:MMw))
  IF (dir == 2) CALL write_hdf_infoREAL(dims,.FALSE.,.FALSE.,'VectorY',array=y2p(SSw:MMw))
  IF (dir == 3) CALL write_hdf_infoREAL(dims,.FALSE.,.FALSE.,'VectorZ',array=y3p(SSw:MMw))
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_1D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  SS1
  INTEGER, INTENT(in)      ::  SS2
  INTEGER, INTENT(in)      ::  SS3
  
  INTEGER, INTENT(in)      ::  NN1
  INTEGER, INTENT(in)      ::  NN2
  INTEGER, INTENT(in)      ::  NN3
  
  INTEGER, INTENT(in)      ::  vel_dir
  
  REAL   , INTENT(out)     ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)        ::  offset_mem(1:3)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  INTEGER                  ::  S3w, M3w, dim3
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Hinweis: phi(1:10,15,12) waere somit                                                                     !
  !       h5dump -d /dsetname -s "11,14,0" -c "1,1,10" filename.h5                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Alt, read2_hdf kann deutlich mehr!                                                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !==========================================================================================================
  CALL filespace_props(vel_dir,(/1,1,1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  dims_mem   = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  dims_data  = (/(NN1-SS1+1   ),(NN2-SS2+1   ),(NN3-SS3+1   )/)
  offset_mem = (/(SS1-b1L     ),(SS2-b2L     ),(SS3-b3L     )/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 lesen ============================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time'             ,scalar=time             )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'dtime'            ,scalar=dtime            )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL read_hdf_infoREAL(1,.TRUE.,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'timestep'         ,scalar=timestep         )
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'write_count'      ,scalar=write_count      )
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL read_hdf_infoINT (1,.TRUE.,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  CALL read_hdf_infoLOG (1,.TRUE.,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  !-----------------------------------------------------------------------------------------------------------
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  CALL h5dclose_f(dset_id  ,herror)
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  ! TEST!!! neu und ungetestet:
  !===========================================================================================================
  !=== Ghost-cell update =====================================================================================
  !===========================================================================================================
  CALL exchange(1,vel_dir,phi)
  CALL exchange(2,vel_dir,phi)
  CALL exchange(3,vel_dir,phi)
  !===========================================================================================================
  
  
  END SUBROUTINE read_hdf
  

  
  !> subroutine to read from an hdf5 file. It is an extended version on routine read_hdf5. 
  !! It can arbitrarily order the read field in space.
  !! @todo: ghost cell update untested.
  SUBROUTINE read2_hdf(filename,dsetname,SS1,SS2,SS3,NN1,NN2,NN3,vel_dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename !< filepath
  CHARACTER(*), INTENT(in) ::  dsetname !< name of the dataset in the hdf5 file
  
  INTEGER, INTENT(in)      ::  SS1      !< lower bounds of dimension 1 
  INTEGER, INTENT(in)      ::  SS2      !< lower bounds of dimension 2
  INTEGER, INTENT(in)      ::  SS3      !< lower bounds of dimension 3
  
  INTEGER, INTENT(in)      ::  NN1      !< upper bounds of dimension 1
  INTEGER, INTENT(in)      ::  NN2      !< upper bounds of dimension 2
  INTEGER, INTENT(in)      ::  NN3      !< upper bounds of dimension 3
  
  INTEGER, INTENT(in)      ::  vel_dir  !< direction of the velocity 0->p, 1,2,3->u,v,w
  
  REAL   , INTENT(out)     ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< requested field
  
  INTEGER(HSIZE_T )        ::  dims_mem(1:3), dims_data(1:3)
  INTEGER(HSSIZE_T)        ::  offset_mem(1:3)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  INTEGER                  ::  S3w, M3w, dim3
  
  INTEGER                  ::  S1r, N1r, i0, iGrid
  INTEGER                  ::  S2r, N2r, j0, jGrid
  INTEGER                  ::  S3r, N3r, k0, kGrid
  
  INTEGER                  ::  Siw(1:3), Miw(1:3)
  
  LOGICAL                ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Hinweis: phi(1:10,15,12) waere somit                                                                     !
  !       h5dump -d /dsetname -s "11,14,0" -c "1,1,10" filename.h5                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Kann im Gegensatz zu "read_hdf" das eingelesene Feld beliebg im Raum anordnen.            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  CALL read_hdf_infoINT (3,.FALSE.,attr_yes,'S1w S2w S3w'      ,array =Siw              )
  CALL read_hdf_infoINT (3,.FALSE.,attr_yes,'M1w M2w M3w'      ,array =Miw              )
  !--- ddemarinis: DA additions
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
 
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  
  !************************************************
  ! Raeumliche Verschiebung des einzulesenden Feldes:
  i0 = 0
  j0 = 0
  k0 = 0
  !************************************************
  
  S1w = Siw(1)
  S2w = Siw(2)
  S3w = Siw(3)
  
  M1w = Miw(1)
  M2w = Miw(2)
  M3w = Miw(3)
  
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  dim3 = M3w-S3w+1
  
  dims_file = (/dim1,dim2,dim3/)
  dims_mem  = (/(N1+b1U-b1L+1),(N2+b2U-b2L+1),(N3+b3U-b3L+1)/)
  
  !-----------------------------------------------------------------------------------------------------------
  iGrid = 0
  jGrid = 0
  kGrid = 0
  IF (vel_dir == 1 .AND. iB(1,1) .GT. 1) THEN
     IF (BC_1L_global  == -2 .AND. ls1 ==  0) iGrid = 1
     IF (BC_1L_global .GT. 0 .AND. ls1 ==  0) iGrid = 2
     IF (BC_1L_global .GT. 0 .AND. ls1 == -1) iGrid = 1
  END IF
  IF (vel_dir == 2 .AND. iB(2,1) .GT. 1) THEN
     IF (BC_2L_global  == -2 .AND. ls2 ==  0) jGrid = 1
     IF (BC_2L_global .GT. 0 .AND. ls2 ==  0) jGrid = 2
     IF (BC_2L_global .GT. 0 .AND. ls2 == -1) jGrid = 1
  END IF
  IF (vel_dir == 3 .AND. iB(3,1) .GT. 1) THEN
     IF (BC_3L_global  == -2 .AND. ls3 ==  0) kGrid = 1
     IF (BC_3L_global .GT. 0 .AND. ls3 ==  0) kGrid = 2
     IF (BC_3L_global .GT. 0 .AND. ls3 == -1) kGrid = 1
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S1w+i0) .LE. (NN1+iShift) .AND. (M1w+i0) .GE. (SS1+iShift)) THEN
     IF ((S1w+i0) .GE. (SS1+iShift)) THEN
        S1r = S1w+i0-iShift
     ELSE
        S1r = SS1
     END IF
     IF ((M1w+i0) .LE. (NN1+iShift)) THEN
        N1r = M1w+i0-iShift
     ELSE
        N1r = NN1
     END IF
     
     dims_data  (1) = N1r-S1r+1
     offset_mem (1) = S1r-b1L
     offset_file(1) = iShift+iGrid-i0
     
     IF (offset_file(1) .LT. 0   ) offset_file(1) = 0
     IF (offset_file(1) .GT. dim1) offset_file(1) = dim1
  ELSE
     dims_data  (1) = 0
     offset_mem (1) = 0
     offset_file(1) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S2w+j0) .LE. (NN2+jShift) .AND. (M2w+j0) .GE. (SS2+jShift)) THEN
     IF ((S2w+j0) .GE. (SS2+jShift)) THEN
        S2r = S2w+j0-jShift
     ELSE
        S2r = SS2
     END IF
     IF ((M2w+j0) .LE. (NN2+jShift)) THEN
        N2r = M2w+j0-jShift
     ELSE
        N2r = NN2
     END IF
     
     dims_data  (2) = N2r-S2r+1
     offset_mem (2) = S2r-b2L
     offset_file(2) = jShift+jGrid-j0
     
     IF (offset_file(2) .LT. 0   ) offset_file(2) = 0
     IF (offset_file(2) .GT. dim2) offset_file(2) = dim2
  ELSE
     dims_data  (2) = 0
     offset_mem (2) = 0
     offset_file(2) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S3w+k0) .LE. (NN3+kShift) .AND. (M3w+k0) .GE. (SS3+kShift)) THEN
     IF ((S3w+k0) .GE. (SS3+kShift)) THEN
        S3r = S3w+k0-kShift
     ELSE
        S3r = SS3
     END IF
     IF ((M3w+k0) .LE. (NN3+kShift)) THEN
        N3r = M3w+k0-kShift
     ELSE
        N3r = NN3
     END IF
     
     dims_data  (3) = N3r-S3r+1
     offset_mem (3) = S3r-b3L
     offset_file(3) = kShift+kGrid-k0
     
     IF (offset_file(3) .LT. 0   ) offset_file(3) = 0
     IF (offset_file(3) .GT. dim3) offset_file(3) = dim3
  ELSE
     dims_data  (3) = 0
     offset_mem (3) = 0
     offset_file(3) = 0
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(3,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  ! TEST!!! neu und ungetestet:
  !===========================================================================================================
  !=== Ghost-cell update =====================================================================================
  !===========================================================================================================
  CALL exchange(1,vel_dir,phi)
  CALL exchange(2,vel_dir,phi)
  CALL exchange(3,vel_dir,phi)
  !===========================================================================================================
  
  
  END SUBROUTINE read2_hdf
  
  
  
  
  
  
  
  
  
  
  ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
  SUBROUTINE read2_2D_hdf(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  N1
  INTEGER, INTENT(in)      ::  N2
  
  INTEGER, INTENT(in)      ::  SS1
  INTEGER, INTENT(in)      ::  SS2
  
  INTEGER, INTENT(in)      ::  NN1
  INTEGER, INTENT(in)      ::  NN2
  
  INTEGER, INTENT(in)      ::  iShift
  INTEGER, INTENT(in)      ::  jShift
  
  INTEGER, INTENT(in)      ::  dir
  
  REAL   , INTENT(out)     ::  phi(1:N1,1:N2)
  
  INTEGER(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
  INTEGER(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  INTEGER                  ::  S1w, M1w, dim1
  INTEGER                  ::  S2w, M2w, dim2
  
  INTEGER                  ::  S1r, N1r, i0
  INTEGER                  ::  S2r, N2r, j0
  
  INTEGER                  ::  Siw(1:2), Miw(1:2)
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen:                                                                                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'dtime'            ,scalar=dtime            )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL read_hdf_infoREAL(1,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'timestep'         ,scalar=timestep         )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL read_hdf_infoINT (1,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  CALL read_hdf_infoLOG (1,.TRUE. ,attr_yes,'new_dtime'        ,scalar=new_dtime        )
  CALL read_hdf_infoINT (2,.FALSE.,attr_yes,'S1w S2w'          ,array =Siw              )
  CALL read_hdf_infoINT (2,.FALSE.,attr_yes,'M1w M2w'          ,array =Miw              )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  
  !************************************************
  ! Raeumliche Verschiebung des einzulesenden Feldes:
  i0 = 0
  j0 = 0
  !************************************************
  
  S1w = Siw(1)
  S2w = Siw(2)
  
  M1w = Miw(1)
  M2w = Miw(2)
  
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  
  dims_file = (/dim1,dim2/)
  dims_mem  = (/N1  ,N2  /)
  
  !-----------------------------------------------------------------------------------------------------------
  IF ((S1w+i0) .LE. (NN1+iShift) .AND. (M1w+i0) .GE. (SS1+iShift)) THEN
     IF ((S1w+i0) .GE. (SS1+iShift)) THEN
        S1r = S1w+i0-iShift
     ELSE
        S1r = SS1
     END IF
     IF ((M1w+i0) .LE. (NN1+iShift)) THEN
        N1r = M1w+i0-iShift
     ELSE
        N1r = NN1
     END IF
     
     dims_data  (1) = N1r-S1r+1
     offset_mem (1) = S1r-1
     offset_file(1) = iShift-i0
     
     IF (offset_file(1) .LT. 0   ) offset_file(1) = 0
     IF (offset_file(1) .GT. dim1) offset_file(1) = dim1
  ELSE
     dims_data  (1) = 0
     offset_mem (1) = 0
     offset_file(1) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF ((S2w+j0) .LE. (NN2+jShift) .AND. (M2w+j0) .GE. (SS2+jShift)) THEN
     IF ((S2w+j0) .GE. (SS2+jShift)) THEN
        S2r = S2w+j0-jShift
     ELSE
        S2r = SS2
     END IF
     IF ((M2w+j0) .LE. (NN2+jShift)) THEN
        N2r = M2w+j0-jShift
     ELSE
        N2r = NN2
     END IF
     
     dims_data  (2) = N2r-S2r+1
     offset_mem (2) = S2r-1
     offset_file(2) = jShift-j0
     
     IF (offset_file(2) .LT. 0   ) offset_file(2) = 0
     IF (offset_file(2) .GT. dim2) offset_file(2) = dim2
  ELSE
     dims_data  (2) = 0
     offset_mem (2) = 0
     offset_file(2) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! ACTHUNG!!! Nicht ganz klar, ob dieses Vorgehen im Sinne des Erfinders ist:
  IF (.NOT. ((dir == -1 .AND. iB(1,1) == 1  ) .OR. (dir == -2 .AND. iB(2,1) == 1  ) .OR. (dir == -3 .AND. iB(3,1) == 1   ) .OR.    &
        &    (dir ==  1 .AND. iB(1,1) == NB1) .OR. (dir ==  2 .AND. iB(2,1) == NB2) .OR. (dir ==  3 .AND. iB(3,1) == NB3))) THEN
     dims_data   = (/0,0/)
     offset_mem  = (/0,0/)
     offset_file = (/0,0/)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(2,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  ! TEST!!! fehlt noch:
  !===========================================================================================================
  !=== Ghost-cell update =====================================================================================
  !===========================================================================================================
  !CALL exchange(1,vel_dir,phi)
  !CALL exchange(2,vel_dir,phi)
  !CALL exchange(3,vel_dir,phi)
  !===========================================================================================================
  
  
  END SUBROUTINE read2_2D_hdf
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read2_1D_hdf(filename,dsetname,NN,SSp,NNp,iShift,dir,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  NN
  INTEGER, INTENT(in)      ::  SSp
  INTEGER, INTENT(in)      ::  NNp
  INTEGER, INTENT(in)      ::  iShift
  INTEGER, INTENT(in)      ::  dir
  
  REAL   , INTENT(out)     ::  phi(1:NN)
  
  INTEGER(HSIZE_T )        ::  dims_file  (1), dims_mem  (1), dims_data(1)
  INTEGER(HSSIZE_T)        ::  offset_file(1), offset_mem(1)
  
  INTEGER                  ::  SSw, MMw, dims
  INTEGER                  ::  SSr, NNr, i0
  
  LOGICAL                  ::  attr_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen:                                                                                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Attribute lesen =======================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open dataset '//dsetname//' !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! Read the attributes:
  attr_yes = .TRUE.
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'SSw'           ,scalar=SSw           ) ! TEST!!! SSw, MMw unschoen ...
  CALL read_hdf_infoINT (1     ,.TRUE. ,attr_yes,'MMw'           ,scalar=MMw           )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5dclose_f(dset_id,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  
  !************************************************
  ! Raeumliche Verschiebung des einzulesenden Feldes:
  i0 = 0
  !************************************************
  
  dims = MMw-SSw+1
  
  dims_file = (/dims/)
  dims_mem  = (/NN /)
  
  !-----------------------------------------------------------------------------------------------------------
  IF ((SSw+i0) .LE. (NNp+iShift) .AND. (MMw+i0) .GE. (SSp+iShift)) THEN
     IF ((SSw+i0) .GE. (SSp+iShift)) THEN
        SSr = SSw+i0-iShift
     ELSE
        SSr = SSp
     END IF
     IF ((MMw+i0) .LE. (NNp+iShift)) THEN
        NNr = MMw+i0-iShift
     ELSE
        NNr = NNp
     END IF
     
     dims_data  (1) = NNr-SSr+1
     offset_mem (1) = SSr-1
     offset_file(1) = iShift-i0
     
     IF (offset_file(1) .LT. 0   ) offset_file(1) = 0
     IF (offset_file(1) .GT. dims) offset_file(1) = dims
  ELSE
     dims_data  (1) = 0
     offset_mem (1) = 0
     offset_file(1) = 0
  END IF
  !-----------------------------------------------------------------------------------------------------------
  ! TEST!!! alle Prozesse lesen!
  !IF (.NOT. ((dir == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) .OR.     &
  !      &    (dir == 2 .AND. iB(1,1) == 1 .AND. iB(3,1) == 1) .OR.     &
  !      &    (dir == 3 .AND. iB(1,1) == 1 .AND. iB(2,1) == 1))) THEN
  !   dims_data   = (/0/)
  !   offset_mem  = (/0/)
  !   offset_file = (/0/)
  !END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(1,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE read2_1D_hdf
  
  
  
  
  
  
  
  
  
  
  

  SUBROUTINE write_stats_hdf_4D(filename,dsetname, SS1,SS2,SS3,SS4, NN1,NN2,NN3,NN4, phi)

  ! ---mjohn 230412

  ! function is suited to write 4D data, specifically designed for spatial mode analysis
  ! dims 1 and 2 are spatial directions, dim1 is x1, dim2 is x2, data is COLLOCATED (count from SSi to NNi).
  ! dims 3 and 4 are integer values, typically Fourier and Hermite modes (count from SSi to NNi).

  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib
  USE HDF5
  USE MPI

  IMPLICIT NONE

  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname

  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3       ! supposed to be 1 (1:NN3)
  INTEGER, INTENT(in)    ::  SS4       ! supposed to be 0 (0:NN4)

  INTEGER, INTENT(in)    ::  NN1       ! supposed to be size in 1 direction (wall-normal)
  INTEGER, INTENT(in)    ::  NN2       ! supposed to be size in 2 direction (spanwise)
  INTEGER, INTENT(in)    ::  NN3       ! supposed to be n_frequencymodes (1:NN3)
  INTEGER, INTENT(in)    ::  NN4       ! supposed to be n_Hermitemodes-1 (0:NN4)

  REAL   , INTENT(in)    ::  phi(SS1:NN1, SS2:NN2, SS3:NN3, SS4:NN4)

  INTEGER                ::  Fourier(SS3:NN3), Hermite(SS4:NN4)
  INTEGER                ::  stride(4)
  INTEGER(HSIZE_T )      ::  dims_file  (4), dims_mem  (4), dims_data(4), stride_mem(4)
  INTEGER(HSSIZE_T)      ::  offset_file(4), offset_mem(4)

  LOGICAL                ::  attr_yes
  INTEGER                ::  i


  ! this function is for stride == 1 only
  stride(1:4) = 1

  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================


  !===========================================================================================================
  !=== Dimensionen und Offsets (hard-coded, da nur eine Struktur von Analyse ausgeschrieben wird) ============
  !===========================================================================================================
  dims_file  = (/ M1,        M2,        NN3-SS3+1, NN4-SS4+1 /)
  dims_mem   = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)
  dims_data  = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)

  offset_mem = (/ SS1-1,     SS2-1,     0,         0         /)
  offset_file= (/ iShift,    jShift,    0,         0         /)
  stride_mem = (/ 1,         1,         1,         1         /)


  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)

  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(4, dims_file, filespace, herror)
  CALL h5screate_simple_f(4, dims_mem , memspace , herror)

  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_count'      ,scalar=write_count      )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_kalm_count' ,scalar=write_kalm_count )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )

  !-----------------------------------------------------------------------------------------------------------
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror,stride_mem)

  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)

  do i = SS3, NN3
     Fourier(i) = i
  end do
  do i = SS4, NN4
     Hermite(i) = i
  end do

  CALL write_hdf_infoREAL(M1,       .FALSE.,.FALSE.,'VectorX',  array=y1p(SS1:NN1)     )
  CALL write_hdf_infoREAL(M2,       .FALSE.,.FALSE.,'VectorY',  array=y2p(SS2:NN2)     )
  CALL write_hdf_infoINT (NN3-SS3+1,.FALSE.,.FALSE.,'VectorDFT',array=Fourier(SS3:NN3) )
  CALL write_hdf_infoINT (NN4-SS4+1,.FALSE.,.FALSE.,'VectorHe', array=Hermite(SS4:NN4) )
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================


  END SUBROUTINE write_stats_hdf_4D




  SUBROUTINE read_stats_hdf_4D(filename,dsetname, SS1,SS2,SS3,SS4, NN1,NN2,NN3,NN4, phi)


  ! ---mjohn 230412

  ! function is suited to read 4D data, specifically designed for spatial mode analysis
  ! dims 1 and 2 are spatial directions, dim1 is x1, dim2 is x2, data is COLLOCATED (count from SSi to NNi).
  ! dims 3 and 4 are integer values, typically Fourier and Hermite modes (count from SSi to NNi).


  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib
  USE HDF5
  USE MPI

  IMPLICIT NONE

  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname

  INTEGER, INTENT(in)    ::  SS1
  INTEGER, INTENT(in)    ::  SS2
  INTEGER, INTENT(in)    ::  SS3       ! supposed to be 1 (1:NN3)
  INTEGER, INTENT(in)    ::  SS4       ! supposed to be 0 (0:NN4)

  INTEGER, INTENT(in)    ::  NN1       ! supposed to be size in 1 direction (wall-normal)
  INTEGER, INTENT(in)    ::  NN2       ! supposed to be size in 2 direction (spanwise)
  INTEGER, INTENT(in)    ::  NN3       ! supposed to be n_frequencymodes (1:NN3)
  INTEGER, INTENT(in)    ::  NN4       ! supposed to be n_Hermitemodes-1 (0:NN4)

  REAL   , INTENT(out)    ::  phi(SS1:NN1, SS2:NN2, SS3:NN3, SS4:NN4)

  INTEGER(HSIZE_T )      ::  dims_file  (4), dims_mem  (4), dims_data(4)
  INTEGER(HSSIZE_T)      ::  offset_file(4), offset_mem(4)



  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)

  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)

  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF

  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================

  !===========================================================================================================
  !=== Dimensionen und Offsets (hard-coded, da nur eine Struktur von Analyse ausgeschrieben wird) ============
  !===========================================================================================================
  dims_file  = (/ M1,          M2,          NN3-SS3+1, NN4-SS4+1 /)
  dims_mem   = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)
  dims_data  = (/ NN1-SS1+1, NN2-SS2+1, NN3-SS3+1, NN4-SS4+1 /)

  offset_mem = (/ SS1-1,       SS2-1,       0,         0         /)
  offset_file= (/ iShift,      jShift,      0,         0         /)


  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)

  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(4,dims_mem,memspace ,herror)

  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)

  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)

  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)

  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================

  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================


  END SUBROUTINE read_stats_hdf_4D


  ! TEST!!! evtl. schoener: write_hdf_2D??
  ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
  SUBROUTINE write_stats_hdf_2D(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  N1
  INTEGER, INTENT(in)      ::  N2
  
  INTEGER, INTENT(in)      ::  SS1
  INTEGER, INTENT(in)      ::  SS2
  
  INTEGER, INTENT(in)      ::  NN1
  INTEGER, INTENT(in)      ::  NN2
  
  INTEGER, INTENT(in)      ::  iShift
  INTEGER, INTENT(in)      ::  jShift
  
  REAL   , INTENT(in)      ::  phi(SS1:NN1,SS2:NN2)
  
  LOGICAL                  ::  attr_yes
  
  INTEGER(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
  INTEGER(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen vel_dir nicht global        !
  !                eingefuehrt.                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Attribut-Datentyp =====================================================================================
  !===========================================================================================================
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE ,memtypeREAL,herror)
  CALL h5tcopy_f(H5T_NATIVE_INTEGER,memtypeINT ,herror)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  dims_file   = (/ N1 , N2 /)
  dims_mem    = (/(NN1-SS1+1),(NN2-SS2+1)/)
  dims_data   = (/(NN1-SS1+1),(NN2-SS2+1)/)
  offset_mem  = (/ SS1-1, SS2-1/)
  offset_file = (/iShift,jShift/)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== HDF5 schreiben ========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access:
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Create the file collectively:
  CALL h5fcreate_f(filename//'.h5',H5F_ACC_TRUNC_F,file_id,herror,access_prp=plist_id)
  CALL h5pclose_f(plist_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Create file space / memory space:
  CALL h5screate_simple_f(2,dims_file,filespace,herror)
  CALL h5screate_simple_f(2,dims_mem ,memspace ,herror)
  
  ! Create the dataset:
  CALL h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,herror)
  !-----------------------------------------------------------------------------------------------------------
  ! Write the attributes:
  attr_yes = .TRUE.
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'restart'          ,scalar=restart          )
  CALL write_hdf_infoINT (1     ,.TRUE. ,attr_yes,'write_stats_count',scalar=write_stats_count)
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_vect'    ,scalar=time_out_vect    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_scal'    ,scalar=time_out_scal    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time_out_kalm'    ,scalar=time_out_kalm    )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_vect'   ,scalar=dtime_out_vect   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_scal'   ,scalar=dtime_out_scal   )
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'dtime_out_kalm'   ,scalar=dtime_out_kalm   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_vect'   ,scalar=write_out_vect   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_scal'   ,scalar=write_out_scal   )
  CALL write_hdf_infoLOG (1     ,.TRUE. ,attr_yes,'write_out_kalm'   ,scalar=write_out_kalm   )
  
  CALL write_hdf_infoINT (2     ,.FALSE.,attr_yes,'SS1 SS2'          ,array =(/SS1,SS2/)      )
  CALL write_hdf_infoINT (2     ,.FALSE.,attr_yes,'NN1 NN2'          ,array =(/NN1,NN2/)      )
  
  CALL write_hdf_infoREAL(1     ,.TRUE. ,attr_yes,'time'             ,scalar=time             )
  !-----------------------------------------------------------------------------------------------------------

  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset write:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Write the dataset collectively:
  CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5dclose_f(dset_id  ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  !-----------------------------------------------------------------------------------------------------------
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  
  END SUBROUTINE write_stats_hdf_2D
  
  ! Anmerkung: Unterscheidet sich im Wesentlichen durch die Dimensionen und den Rang von phi ...
  SUBROUTINE read_stats_hdf_2D(filename,dsetname,N1,N2,SS1,SS2,NN1,NN2,iShift,jShift,phi)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(in) ::  filename
  CHARACTER(*), INTENT(in) ::  dsetname
  
  INTEGER, INTENT(in)      ::  N1
  INTEGER, INTENT(in)      ::  N2
  
  INTEGER, INTENT(in)      ::  SS1
  INTEGER, INTENT(in)      ::  SS2
  
  INTEGER, INTENT(in)      ::  NN1
  INTEGER, INTENT(in)      ::  NN2
  
  INTEGER, INTENT(in)      ::  iShift
  INTEGER, INTENT(in)      ::  jShift
  
  REAL   , INTENT(out)     ::  phi(SS1:NN1,SS2:NN2)
  
  INTEGER(HSIZE_T )        ::  dims_file  (1:2), dims_mem  (1:2), dims_data(1:2)
  INTEGER(HSSIZE_T)        ::  offset_file(1:2), offset_mem(1:2)
  
  !===========================================================================================================
  !=== File oeffnen ==========================================================================================
  !===========================================================================================================
  ! Setup file access property list with parallel I/O access
  CALL h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,herror)
  CALL h5pset_fapl_mpio_f(plist_id,COMM_CART,MPI_INFO_NULL,herror)
  
  ! Open the file collectively
  CALL h5fopen_f(filename//'.h5',H5F_ACC_RDWR_F,file_id,herror,access_prp=plist_id)
  
  IF (herror == -1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open file '//filename//'.h5 !'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  CALL h5pclose_f(plist_id,herror)
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Dimensionen und Offsets ===============================================================================
  !===========================================================================================================
  dims_file   = (/N1  ,N2  /)
  dims_mem    = (/(NN1-SS1+1),(NN2-SS2+1)/)
  dims_data   = (/(NN1-SS1+1),(NN2-SS2+1)/)
  offset_mem  = (/ SS1-1, SS2-1/)
  offset_file = (/iShift,jShift/)
  !===========================================================================================================
  
  !===========================================================================================================
  !=== Feld lesen ============================================================================================
  !===========================================================================================================
  phi = 0.
  !-----------------------------------------------------------------------------------------------------------
  ! Open the dataset:
  CALL h5dopen_f(file_id,dsetname,dset_id,herror)
  
  ! Get file space / create memory space:
  CALL h5dget_space_f    (dset_id   ,filespace,herror)
  CALL h5screate_simple_f(2,dims_mem,memspace ,herror)
  
  ! Select hyperslab in the file space / memory space:
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset_file,dims_data,herror)
  CALL h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,offset_mem ,dims_data,herror)
  
  ! Create property list for collective dataset read:
  CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,herror)
  CALL h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,herror)
  
  ! Read the dataset collectively:
  CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,phi,dims_file,herror,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
  
  CALL h5pclose_f(plist_id ,herror)
  CALL h5sclose_f(memspace ,herror)
  CALL h5sclose_f(filespace,herror)
  
  CALL h5dclose_f(dset_id  ,herror)
  !===========================================================================================================
  
  !===========================================================================================================
  !=== File schliessen =======================================================================================
  !===========================================================================================================
  CALL h5fclose_f(file_id,herror)
  !===========================================================================================================
  
  END SUBROUTINE read_stats_hdf_2D
  


  SUBROUTINE write_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  NN
  LOGICAL          , INTENT(in)  ::  scalar_yes
  LOGICAL          , INTENT(in)  ::  attr_yes
  CHARACTER(*)     , INTENT(in)  ::  name
  REAL   , OPTIONAL, INTENT(in)  ::  array(1:NN)
  REAL   , OPTIONAL, INTENT(in)  ::  scalar
  
  INTEGER(HID_T)                 ::  memspace
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  REAL                           ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (scalar_yes) THEN
     value = scalar
  ELSE
     value = array
  END IF
  
  CALL h5screate_simple_f(1,dim_mem,memspace,herror)
  
  IF (attr_yes) THEN
     CALL h5acreate_f(dset_id,name,memtypeREAL,memspace,attr_id,herror)
     CALL h5awrite_f(attr_id,memtypeREAL,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,memspace,attr_id,herror)
     CALL H5dwrite_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror,memspace)!,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id ,herror)
  END IF
  
  CALL h5sclose_f(memspace,herror)
  
  
  END SUBROUTINE write_hdf_infoREAL
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  NN
  LOGICAL          , INTENT(in)  ::  scalar_yes
  LOGICAL          , INTENT(in)  ::  attr_yes
  CHARACTER(*)     , INTENT(in)  ::  name
  INTEGER, OPTIONAL, INTENT(in)  ::  array(1:NN)
  INTEGER, OPTIONAL, INTENT(in)  ::  scalar
  
  INTEGER(HID_T)                 ::  memspace
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (scalar_yes) THEN
     value = scalar
  ELSE
     value = array
  END IF
  
  CALL h5screate_simple_f(1,dim_mem,memspace,herror)
  
  IF (attr_yes) THEN
     CALL h5acreate_f(dset_id,name,memtypeINT ,memspace,attr_id,herror)
     CALL h5awrite_f(attr_id,memtypeINT ,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,memspace,attr_id,herror)
     CALL H5dwrite_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,memspace)!,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id ,herror)
  END IF
  
  CALL h5sclose_f(memspace,herror)
  
  
  END SUBROUTINE write_hdf_infoINT
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  NN
  LOGICAL          , INTENT(in)  ::  scalar_yes
  LOGICAL          , INTENT(in)  ::  attr_yes
  CHARACTER(*)     , INTENT(in)  ::  name
  LOGICAL, OPTIONAL, INTENT(in)  ::  array(1:NN)
  LOGICAL, OPTIONAL, INTENT(in)  ::  scalar
  
  INTEGER(HID_T)                 ::  memspace
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  m
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (scalar_yes) THEN
     IF (scalar) THEN
        value = 1
     ELSE
        value = 0
     END IF
  ELSE
     DO m = 1, NN
        IF (array(m)) THEN
           value(m) = 1
        ELSE
           value(m) = 0
        END IF
     END DO
  END IF
  
  CALL h5screate_simple_f(1,dim_mem,memspace,herror)
  
  IF (attr_yes) THEN
     CALL h5acreate_f(dset_id,name,memtypeINT ,memspace,attr_id,herror)
     CALL h5awrite_f(attr_id,memtypeINT ,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,memspace,attr_id,herror)
     CALL H5dwrite_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror,memspace)!,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id ,herror)
  END IF
  
  CALL h5sclose_f(memspace,herror)
  
  
  END SUBROUTINE write_hdf_infoLOG
  
  
  
  
  
  
  
  
  
  
  !> subroutine for recovering real valued attribute data from hdf5 files 
  SUBROUTINE read_hdf_infoREAL(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  NN            !< dimension of the output (1=scalar)
  LOGICAL          , INTENT(in)  ::  scalar_yes    !< whether output is scalar
  LOGICAL          , INTENT(in)  ::  attr_yes      !< ?
  CHARACTER(*)     , INTENT(in)  ::  name          !< attribute name
  REAL   , OPTIONAL, INTENT(out) ::  array(1:NN)   !< multidimensional output
  REAL   , OPTIONAL, INTENT(out) ::  scalar        !< onedimensional output
  
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  REAL                           ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (attr_yes) THEN
     CALL h5aopen_name_f(dset_id,name,attr_id,herror)
     IF (rank == 0) CALL h5aread_f(attr_id,memtypeREAL,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dopen_f(file_id,name,attr_id,herror)
     IF (rank == 0) CALL H5dread_f(attr_id,H5T_NATIVE_DOUBLE,value,dim_mem,herror) !,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id,herror)
  END IF
  
  CALL MPI_BCAST(value,NN,MPI_REAL8,0,COMM_CART,merror)
  
  IF (scalar_yes) THEN
     scalar = value(1)
  ELSE
     array  = value
  END IF
  
  
  END SUBROUTINE read_hdf_infoREAL
  
  
  
  
  
  
  
  
  
  
  !> subroutine for recovering integer valued attribute data from hdf5 files.
  SUBROUTINE read_hdf_infoINT(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  NN         !< dimension of the output
  LOGICAL          , INTENT(in)  ::  scalar_yes !< whether output is scalar
  LOGICAL          , INTENT(in)  ::  attr_yes   !< ?
  CHARACTER(*)     , INTENT(in)  ::  name       !< attribute name
  INTEGER, OPTIONAL, INTENT(out) ::  array(1:NN)!< multidimensional output
  INTEGER, OPTIONAL, INTENT(out) ::  scalar     !< onedimensional output
  
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (attr_yes) THEN
     CALL h5aopen_name_f(dset_id,name,attr_id,herror)
     IF (rank == 0) CALL h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dopen_f(file_id,name,attr_id,herror)
     IF (rank == 0) CALL H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror)!,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id,herror)
  END IF
  
  CALL MPI_BCAST(value,NN,MPI_INTEGER,0,COMM_CART,merror)
  
  IF (scalar_yes) THEN
     scalar = value(1)
  ELSE
     array  = value
  END IF
  
  
  END SUBROUTINE read_hdf_infoINT
  
  
  
  
  
  
  
  
  
  
  !< subroutine for recovering logical valued attribute data from hdf5 files 
  SUBROUTINE read_hdf_infoLOG(NN,scalar_yes,attr_yes,name,array,scalar)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  NN          !< dimension of output
  LOGICAL          , INTENT(in)  ::  scalar_yes  !< whether output is scalar
  LOGICAL          , INTENT(in)  ::  attr_yes    !< ?
  CHARACTER(*)     , INTENT(in)  ::  name        !< attribute name
  LOGICAL, OPTIONAL, INTENT(out) ::  array(1:NN) !< multidimensional output
  LOGICAL, OPTIONAL, INTENT(out) ::  scalar      !< onedimensional output
  
  INTEGER(HID_T)                 ::  attr_id
  INTEGER(HSIZE_T)               ::  dim_mem(1)
  INTEGER                        ::  m
  INTEGER                        ::  value(1:NN)
  
  
  dim_mem = (/NN/)
  
  IF (attr_yes) THEN
     CALL h5aopen_name_f(dset_id,name,attr_id,herror)
     IF (rank == 0) CALL h5aread_f(attr_id,memtypeINT,value,dim_mem,herror)
     CALL h5aclose_f(attr_id,herror)
  ELSE
     CALL h5dopen_f(file_id,name,attr_id,herror)
     IF (rank == 0) CALL H5dread_f(attr_id,H5T_NATIVE_INTEGER,value,dim_mem,herror) !,H5S_ALL_F,H5S_ALL_F,H5P_DEFAULT_F)
     CALL h5dclose_f(attr_id,herror)
  END IF
  
  CALL MPI_BCAST(value,NN,MPI_LOGICAL,0,COMM_CART,merror)
  
  IF (scalar_yes) THEN
     IF (value(1) == 1) THEN
        scalar = .TRUE.
     ELSE
        scalar = .FALSE.
     END IF
  ELSE
     DO m = 1, NN
        IF (value(m) == 1) THEN
           array(m) = .TRUE.
        ELSE
           array(m) = .FALSE.
        END IF
     END DO
  END IF
  
  
  END SUBROUTINE read_hdf_infoLOG
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE filespace_props(vel_dir,stride,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
  
  IMPLICIT NONE
  
  INTEGER          , INTENT(in)  ::  vel_dir
  INTEGER          , INTENT(in)  ::  stride(1:3)
  INTEGER          , INTENT(out) ::  S1w, M1w, dim1
  INTEGER          , INTENT(out) ::  S2w, M2w, dim2
  INTEGER          , INTENT(out) ::  S3w, M3w, dim3
  
  
  !-----------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Intervallgrenzen, Offsets und Block-Groessen werden wegen "vel_dir" nicht global eingefuehrt, so   !
  !                dass diese Routine immer wieder von Neuem aufgerufen werden muss.                                      !
  !              - Routine dient in erster Linie der kuerzeren Schreibweise / besseren Uebersicht                         !
  !              - Default bei Intervallgrenzen sind periodische / Nachbarblock-RB                                        !
  !-----------------------------------------------------------------------------------------------------------------------!
  
  
  S1w = 2 + ls1
  S2w = 2 + ls2
  S3w = 2 + ls3
  
  M1w = (M1-1)/stride(1) + 1 + ls1
  M2w = (M2-1)/stride(2) + 1 + ls2
  M3w = (M3-1)/stride(3) + 1 + ls3
  
  offset_file = (/iShift/stride(1),jShift/stride(2),kShift/stride(3)/)
  
  
  !===========================================================================================================
  IF (BC_1L_global /= -1) THEN
     IF (vel_dir == 1) THEN
        IF (BC_1L_global == -2) THEN
           S1w = 1
           IF (iB(1,1) .GT. 1 .AND. ls1 ==  0) offset_file(1) = offset_file(1) + 1
        ELSE
           S1w = 0
           IF (iB(1,1) .GT. 1 .AND. ls1 ==  0) offset_file(1) = offset_file(1) + 2
           IF (iB(1,1) .GT. 1 .AND. ls1 == -1) offset_file(1) = offset_file(1) + 1
        END IF
        IF (BC_1U_global == -2) THEN
           M1w = M1-1
        ELSE
           M1w = M1
        END IF
     ELSE
        S1w = 1
        M1w = (M1-1)/stride(1) + 1
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L_global /= -1) THEN
     IF (vel_dir == 2) THEN
        IF (BC_2L_global == -2) THEN
           S2w = 1
           IF (iB(2,1) .GT. 1 .AND. ls2 ==  0) offset_file(2) = offset_file(2) + 1
        ELSE
           S2w = 0
           IF (iB(2,1) .GT. 1 .AND. ls2 ==  0) offset_file(2) = offset_file(2) + 2
           IF (iB(2,1) .GT. 1 .AND. ls2 == -1) offset_file(2) = offset_file(2) + 1
        END IF
        IF (BC_2U_global == -2) THEN
           M2w = M2-1
        ELSE
           M2w = M2
        END IF
     ELSE
        S2w = 1
        M2w = (M2-1)/stride(2) + 1
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_3L_global /= -1) THEN
     IF (vel_dir == 3) THEN
        IF (BC_3L_global == -2) THEN
           S3w = 1
           IF (iB(3,1) .GT. 1 .AND. ls3 ==  0) offset_file(3) = offset_file(3) + 1
        ELSE
           S3w = 0
           IF (iB(3,1) .GT. 1 .AND. ls3 ==  0) offset_file(3) = offset_file(3) + 2
           IF (iB(3,1) .GT. 1 .AND. ls3 == -1) offset_file(3) = offset_file(3) + 1
        END IF
        IF (BC_3U_global == -2) THEN
           M3w = M3-1
        ELSE
           M3w = M3
        END IF
     ELSE
        S3w = 1
        M3w = (M3-1)/stride(3) + 1
     END IF
  END IF
  !===========================================================================================================
  
  dim1 = M1w-S1w+1
  dim2 = M2w-S2w+1
  dim3 = M3w-S3w+1
  
  dims_file = (/dim1,dim2,dim3/)
  
  
  END SUBROUTINE filespace_props
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE lambda
  
  IMPLICIT NONE
  
  REAL                   ::  Dvel(1:3,1:3)
  REAL                   ::  TT(1:6)
  
  REAL                   ::  PP, QQ, RR
  REAL                   ::  rho, theta
  REAL                   ::  eps
  
  REAL                   ::  temp
  REAL                   ::  lam(1:3)
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m, n
  REAL                   ::  pi
  
  
  !-----------------------------------------------------------------!
  ! Anmerkung: [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]    !
  !-----------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  eps = 10.**(-5)
  
  !===========================================================================================================
  !=== Geschwindigkeiten interpolieren =======================================================================
  !===========================================================================================================
  
  ! TEST!!! kann man auch weglassen (wird schon vorher ausgefuehrt), bzw. exchange-Parameter einfuehren!
  CALL exchange(1,1,vel(b1L,b2L,b3L,1))
  CALL exchange(2,2,vel(b1L,b2L,b3L,2))
  CALL exchange(3,3,vel(b1L,b2L,b3L,3))
  
  !===========================================================================================================
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
  END DO
  
  !===========================================================================================================
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
  END DO
  
  !===========================================================================================================
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,3) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              nl(i,j,k,3) = nl(i,j,k,3) + cIwp(kk,k)*vel(i,j,k+kk,3)
           END DO
        END DO
     END DO
  END DO
  
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== lambda ================================================================================================
  !===========================================================================================================
  
  CALL exchange_all_all(.FALSE.,nl)
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           
           Dvel = 0.
           
           !--- d/dx -----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO ii = b1L, b1U
              Dvel(:,1) = Dvel(:,1) + cp1(ii,i)*nl(i+ii,j,k,:)
           END DO
           
           !--- d/dy -----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO jj = b2L, b2U
              Dvel(:,2) = Dvel(:,2) + cp2(jj,j)*nl(i,j+jj,k,:)
           END DO
           
           !--- d/dz -----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO kk = b3L, b3U
              Dvel(:,3) = Dvel(:,3) + cp3(kk,k)*nl(i,j,k+kk,:)
           END DO
           
           
           !--- Tensor (symmetrisch) -------------------------------------------------------------------------
           TT(1) = 2.*(Dvel(1,1)*Dvel(1,1) + Dvel(1,2)*Dvel(2,1) + Dvel(1,3)*Dvel(3,1))
           TT(2) = 2.*(Dvel(2,1)*Dvel(1,2) + Dvel(2,2)*Dvel(2,2) + Dvel(2,3)*Dvel(3,2))
           TT(3) = 2.*(Dvel(3,1)*Dvel(1,3) + Dvel(3,2)*Dvel(2,3) + Dvel(3,3)*Dvel(3,3))
           
           TT(4) = (Dvel(1,1) + Dvel(2,2))*(Dvel(1,2) + Dvel(2,1)) + Dvel(1,3)*Dvel(3,2) + Dvel(2,3)*Dvel(3,1)
           TT(5) = (Dvel(1,1) + Dvel(3,3))*(Dvel(1,3) + Dvel(3,1)) + Dvel(1,2)*Dvel(2,3) + Dvel(3,2)*Dvel(2,1)
           TT(6) = (Dvel(2,2) + Dvel(3,3))*(Dvel(2,3) + Dvel(3,2)) + Dvel(2,1)*Dvel(1,3) + Dvel(3,1)*Dvel(1,2)
           
           
           !--- Invarianten ----------------------------------------------------------------------------------
           PP = TT(1) + TT(2) + TT(3)
           QQ = TT(1)*TT(2) + TT(1)*TT(3) + TT(2)*TT(3) - TT(4)**2 - TT(5)**2 - TT(6)**2
           RR = TT(1)*TT(2)*TT(3) + 2.*TT(4)*TT(5)*TT(6) - TT(1)*TT(6)**2 - TT(2)*TT(5)**2 - TT(3)*TT(4)**2
           
           
           !--- Eigenwerte -----------------------------------------------------------------------------------
           rho = (PP/3.)**2 - QQ/3.
           IF (rho .LE. eps) THEN
              !----------------------------------------------------------------------------------------------!
              ! y:=lam-PP/3.                                                                                 !
              ! y**3+p*y+q=0.                                                                                !
              ! p:=QQ-PP**2/3.                                                                               !
              ! q:=-2*(PP/3)**3+PP*QQ/3.-RR                                                                  !
              ! TT ist symmetrisch ==> lam(1:3) reell <==> D=(p/3.)**3+(q/2.)**2 .LE. 0.                     !
              !                    ==> falls p=QQ-PP**2/3.=0.                                                !
              !                    ==> q:=-2*(PP/3)**3+PP**3/9.-RR = (PP/3)**3-RR                            !
              !                    ==> D=(q/2.)**2 .LE. 0. <==> q=0. <==> RR=(PP/3)**3 <==> lam(1:3)=lam(3)  !
              !                    ==> y=0.                                                                  !
              !                    ==> lam=PP/3.=RR**(1./3.)                                                 !
              !----------------------------------------------------------------------------------------------!
              res(i,j,k) = PP/3.
           ELSE
              rho = SQRT(rho)
              QQ  = ((PP/3.)**3 - QQ*(PP/6.) + RR/2.)/rho**3
              
              IF (ABS(QQ) .LT. 1.) THEN
                 theta = ACOS(QQ)/3.
              ELSE
                 IF (QQ .GT. 0.) THEN
                    theta = 0.
                 ELSE
                    theta = pi/3.
                 END IF
              END IF
              
              lam(1) = COS(theta           )
              lam(2) = COS(theta + 2.*pi/3.)
              lam(3) = COS(theta + 4.*pi/3.)
              
              
              !--- sortieren ---------------------------------------------------------------------------------
              DO m = 1, 3
!pgi$ unroll = n:8
                 DO n = m+1, 3
                    IF (lam(m) .GT. lam(n)) THEN
                       temp   = lam(n)
                       lam(n) = lam(m)
                       lam(m) = temp
                    END IF
                 END DO
              END DO
              
              ! Faktor 1/2, da bei TT(1:6) mit 2 durchmultipliziert wurde ...
              res(i,j,k) = rho*lam(2) + PP/6.
           END IF
           
        END DO
     END DO
  END DO
  
  
  END SUBROUTINE lambda
  

  



  SUBROUTINE vorticity_2D

  IMPLICIT NONE
  
  REAL                   ::  Dvel(1:2)
  
  
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  
  !-----------------------------------------------------------------!
  ! Anmerkung: [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]    !
  !-----------------------------------------------------------------!
  
  !===========================================================================================================
  !=== Geschwindigkeiten interpolieren =======================================================================
  !===========================================================================================================
  
  ! TEST!!! kann man auch weglassen (wird schon vorher ausgefuehrt), bzw. exchange-Parameter einfuehren!
  CALL exchange(1,1,vel(b1L,b2L,b3L,1))
  CALL exchange(2,2,vel(b1L,b2L,b3L,2))
  
  !===========================================================================================================
  
  k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nl(i,j,k,1) = nl(i,j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           END DO
        END DO
     END DO
  
  !===========================================================================================================
  
  k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           nl(i,j,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nl(i,j,k,2) = nl(i,j,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           END DO
        END DO
     END DO
  
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== vorticity =============================================================================================
  !===========================================================================================================
  
  CALL exchange_all_all(.FALSE.,nl)
  
  k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           
           Dvel = 0.
           
           !--- dv/dx ----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO ii = b1L, b1U
              Dvel(1) = Dvel(1) + cp1(ii,i)*nl(i+ii,j,k,2)
           END DO
           
           !--- du/dy ----------------------------------------------------------------------------------------
!pgi$ unroll = n:8
           DO jj = b2L, b2U
              Dvel(2) = Dvel(2) + cp2(jj,j)*nl(i,j+jj,k,1)
           END DO

           res(i,j,k) = Dvel(1) - Dvel(2)

        END DO
     END DO

  END SUBROUTINE vorticity_2D


  SUBROUTINE write_xdmf_xml()

    IMPLICIT NONE

    CHARACTER(LEN=8)     ::  count_char
    INTEGER              ::  S1w, S2w, S3w
    INTEGER              ::  M1w, M2w, M3w
    INTEGER              ::  dim1, dim2, dim3

    IF (rank .EQ. 0) THEN

      !=========================================================================================================
      !=== Write count string ==================================================================================
      !=========================================================================================================
      CALL num_to_string(8,write_count,count_char)
      !=========================================================================================================

      WRITE(*,'(2a)') 'writing XMF file(s) for count_char ...', count_char

      !=========================================================================================================
      !=== Dimensions of fields ================================================================================
      !=========================================================================================================
      ! The way IMPACT currently writes fields the vel_dir argument is always
      ! equal to 0.
      IF (dimens .EQ. 3) THEN
        IF (write_large) THEN
          CALL filespace_props(0,stride_large,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
          CALL write_3D_xdmf('IMPACT_xdmf3d','large',count_char,dim1,dim2,dim3)
        END IF
        IF (write_med  ) THEN
          CALL filespace_props(0,stride_med  ,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
          CALL write_3D_xdmf('IMPACT_xdmf3d','med'  ,count_char,dim1,dim2,dim3)
        END IF
        IF (write_small) THEN
          CALL filespace_props(0,stride_small,S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
          CALL write_3D_xdmf('IMPACT_xdmf3d','small',count_char,dim1,dim2,dim3)
        END IF
      ELSE
        ! will this work for 2D?
        CALL filespace_props(0,(/1, 1, 1/),S1w,S2w,S3w,M1w,M2w,M3w,dim1,dim2,dim3)
        CALL write_2D_xdmf('IMPACT_xdmf2d',count_char,dim1,dim2)
      END IF


    END IF

  END SUBROUTINE write_xdmf_xml


  SUBROUTINE write_3D_xdmf(filebase, stride_string, count_char, dim1, dim2, dim3)

    IMPLICIT NONE

    CHARACTER(*)    , INTENT(in)  ::  filebase, stride_string
    CHARACTER(LEN=8), INTENT(in)  ::  count_char
    CHARACTER(LEN=8)              ::  count_char2

    INTEGER, INTENT(in)           ::  dim1, dim2, dim3
    INTEGER                       ::  xmf=995

    CHARACTER(*), PARAMETER       ::  fmt1 = "(a)"
    CHARACTER(*), PARAMETER       ::  fmt2 = "(a,i0,a)"
    CHARACTER(*), PARAMETER       ::  fmt3 = "(a,i0,1x,i0,1x,i0,a)"
    CHARACTER(*), PARAMETER       ::  fmt4 = "(a,f15.9,a)"
    CHARACTER(*), PARAMETER       ::  fmt5 = "(a,f15.6,a,i0,a)"
    CHARACTER(*), PARAMETER       ::  fmt6 = "(a,f15.6,a,i0,1x,i0,1x,i0,a)"
    CHARACTER(LEN=2)              ::  id,phs
    CHARACTER(LEN=50)             ::  write_dir

    !=== Open file to write into ==============================================================================
    OPEN(xmf, FILE=filebase//'_'//stride_string//'_'//count_char//'.xmf')

    WRITE(xmf,fmt1) '<?xml version="1.0" ?>'
    WRITE(xmf,fmt1) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    WRITE(xmf,fmt1) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    WRITE(xmf,fmt1) '  <Domain>'
    !=== Create grid that holds all data ======================================================================
    WRITE(xmf,fmt1) '    <Grid Name="IMPACT mesh" GridType="Uniform">'
    !=== Define grid's coordinates (rectilinear) ==============================================================
    WRITE(xmf,fmt3) '      <Topology TopologyType="3DRectMesh" Dimensions="',dim3,dim2,dim1,'"/>'
    WRITE(xmf,fmt1) '        <Geometry GeometryType="VXVYVZ">'
    IF (scale_output_yes) WRITE(xmf,fmt5) '          <DataItem ItemType="Function" Function="',L_ref,' * $0" Dimensions="',dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt2) '          <DataItem Format="HDF" Dimensions="',dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//stride_string//'_'//count_char//'.h5:/VectorX'
    WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(xmf,fmt5) '          <DataItem ItemType="Function" Function="',L_ref,' * $0" Dimensions="',dim2,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt2) '          <DataItem Format="HDF" Dimensions="',dim2,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//stride_string//'_'//count_char//'.h5:/VectorY'
    WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(xmf,fmt5) '          <DataItem ItemType="Function" Function="',L_ref,' * $0" Dimensions="',dim3,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt2) '          <DataItem Format="HDF" Dimensions="',dim3,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//stride_string//'_'//count_char//'.h5:/VectorZ'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Geometry>'
    !=== Add information about time ============================================================================
    IF (scale_output_yes) THEN
      WRITE(xmf,fmt4) '        <Time TimeType="Single" Value="',time * L_ref/U_ref,'"/>'
    ELSE
      WRITE(xmf,fmt4) '        <Time TimeType="Single" Value="',time,'"/>'
    END IF
    !=== Add the datasets ======================================================================================
    !--- Pressure
    WRITE(xmf,fmt1) '        <Attribute Name="Pressure" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',rho_fluid * U_ref**2.,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//stride_string//'_'//count_char//'.h5:/pre'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- X-Velocity
    WRITE(xmf,fmt1) '        <Attribute Name="Velocity_X" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            velX_'//stride_string//'_'//count_char//'.h5:/velX'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- Y-Velocity
    WRITE(xmf,fmt1) '        <Attribute Name="Velocity_Y" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            velY_'//stride_string//'_'//count_char//'.h5:/velY'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- Z-Velocity
    WRITE(xmf,fmt1) '        <Attribute Name="Velocity_Z" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            velZ_'//stride_string//'_'//count_char//'.h5:/velZ'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- Vorticity (2D) or Lambda-2 (3D)
    IF (write_lambda2_yes) THEN
      WRITE(xmf,fmt1) '        <Attribute Name="Lambda_2" Center="Node">'
      IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref/L_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            lamb_'//stride_string//'_'//count_char//'.h5:/lamb'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
    END IF
    !--- Forces
    IF (write_force_yes) THEN
      !--- X-Force
      WRITE(xmf,fmt1) '        <Attribute Name="Volume_Force_X" Center="Node">'
      IF (.NOT. scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',L_ref / (rho_fluid * U_ref**2.),' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            forceX_'//stride_string//'_'//count_char//'.h5:/forceX'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (.NOT. scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
      !--- Y-Force
      WRITE(xmf,fmt1) '        <Attribute Name="Volume_Force_Y" Center="Node">'
      IF (.NOT. scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',L_ref / (rho_fluid * U_ref**2.),' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            forceY_'//stride_string//'_'//count_char//'.h5:/forceY'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (.NOT. scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
      !--- Z-Force
      WRITE(xmf,fmt1) '        <Attribute Name="Volume_Force_Z" Center="Node">'
      IF (.NOT. scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',L_ref / (rho_fluid * U_ref**2.),' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            forceZ_'//stride_string//'_'//count_char//'.h5:/forceZ'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (.NOT. scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
    END IF

    !---- Covariance
    IF (write_covariance_yes) THEN
      if (dtime_out_kalm.ne.0.0) then
         phase = mod(write_kalm_count-1,intervals) + 1
         CALL num_to_string(2,phase,phs)
         write_dir = './kf_result/phase_'//phs//'/'
         CALL num_to_string(8,write_kalm_count,count_char2)
         !--- Kalman gain ------------------------------------------------------------------------------------------------------------------------------
         !--- data_XX
         WRITE(xmf,fmt1) '        <Attribute Name="GainX" Center="Node">'
         IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
         WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
         WRITE(xmf,fmt1) '            '//trim(write_dir)//'gainX_phase'//phs//'_'//count_char2//'.h5:/gainX'
         WRITE(xmf,fmt1) '          </DataItem>'
         IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
         WRITE(xmf,fmt1) '        </Attribute>'
         !--- data_YY
         WRITE(xmf,fmt1) '        <Attribute Name="GainY" Center="Node">'
         IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
         WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
         WRITE(xmf,fmt1) '            '//trim(write_dir)//'gainY_phase'//phs//'_'//count_char2//'.h5:/gainY'
         WRITE(xmf,fmt1) '          </DataItem>'
         IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
         WRITE(xmf,fmt1) '        </Attribute>'
         !--- data_ZZ
         WRITE(xmf,fmt1) '        <Attribute Name="GainZ" Center="Node">'
         IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
         WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim3,dim2,dim1,'" NumberType="Float" Precision="8">'
         WRITE(xmf,fmt1) '            '//trim(write_dir)//'gainZ_phase'//phs//'_'//count_char2//'.h5:/gainZ'
         WRITE(xmf,fmt1) '          </DataItem>'
         IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
         WRITE(xmf,fmt1) '        </Attribute>'
      end if
    END IF

    WRITE(xmf,fmt1) '    </Grid>'
    WRITE(xmf,fmt1) '  </Domain>'
    WRITE(xmf,fmt1) '</Xdmf>'

    CLOSE(xmf)

  END SUBROUTINE write_3D_xdmf



  SUBROUTINE write_2D_xdmf(filebase, count_char, dim1, dim2)

    IMPLICIT NONE

    CHARACTER(*)    , INTENT(in)  ::  filebase
    CHARACTER(LEN=8), INTENT(in)  ::  count_char

    INTEGER, INTENT(in)           ::  dim1, dim2
    INTEGER                       ::  xmf=995

    CHARACTER(*), PARAMETER       ::  fmt1 = "(a)"
    CHARACTER(*), PARAMETER       ::  fmt2 = "(a,i0,a)"
    CHARACTER(*), PARAMETER       ::  fmt3 = "(a,i0,1x,i0,a)"
    CHARACTER(*), PARAMETER       ::  fmt4 = "(a,f15.9,a)"
    CHARACTER(*), PARAMETER       ::  fmt5 = "(a,f15.6,a,i0,a)"
    CHARACTER(*), PARAMETER       ::  fmt6 = "(a,f15.6,a,i0,1x,i0,a)"

    !=== Open file to write into ==============================================================================
    OPEN(xmf, FILE=filebase//'_'//count_char//'.xmf')

    WRITE(xmf,fmt1) '<?xml version="1.0" ?>'
    WRITE(xmf,fmt1) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    WRITE(xmf,fmt1) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    WRITE(xmf,fmt1) '  <Domain>'
    !=== Create grid that holds all data ======================================================================
    WRITE(xmf,fmt1) '    <Grid Name="IMPACT mesh" GridType="Uniform">'
    !=== Define grid's coordinates (rectilinear) ==============================================================
    WRITE(xmf,fmt3) '      <Topology TopologyType="2DRectMesh" Dimensions="',dim2,dim1,'"/>'
    WRITE(xmf,fmt1) '        <Geometry GeometryType="VXVY">'
    IF (scale_output_yes) WRITE(xmf,fmt5) '          <DataItem ItemType="Function" Function="',L_ref,' * $0" Dimensions="',dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt2) '          <DataItem Format="HDF" Dimensions="',dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//count_char//'.h5:/VectorX'
    WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(xmf,fmt5) '          <DataItem ItemType="Function" Function="',L_ref,' * $0" Dimensions="',dim2,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt2) '          <DataItem Format="HDF" Dimensions="',dim2,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//count_char//'.h5:/VectorY'
    WRITE(XMF,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Geometry>'
    !=== Add information about time ============================================================================
    IF (scale_output_yes) THEN
      WRITE(xmf,fmt4) '        <Time TimeType="Single" Value="',time * L_ref/U_ref,'"/>'
    ELSE
      WRITE(xmf,fmt4) '        <Time TimeType="Single" Value="',time,'"/>'
    END IF
    !=== Add the datasets ======================================================================================
    !--- Pressure
    WRITE(xmf,fmt1) '        <Attribute Name="Pressure" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',rho_fluid * U_ref**2.,' * $0" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            pre_'//count_char//'.h5:/pre'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- X-Velocity
    WRITE(xmf,fmt1) '        <Attribute Name="Velocity_X" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            velX_'//count_char//'.h5:/velX'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- Y-Velocity
    WRITE(xmf,fmt1) '        <Attribute Name="Velocity_Y" Center="Node">'
    IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref,' * $0" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
    WRITE(xmf,fmt1) '            velY_'//count_char//'.h5:/velY'
    WRITE(xmf,fmt1) '          </DataItem>'
    IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
    WRITE(xmf,fmt1) '        </Attribute>'
    !--- Vorticity (2D) or Lambda-2 (3D)
    IF (write_lambda2_yes) THEN
      WRITE(xmf,fmt1) '        <Attribute Name="Vorticity" Center="Node">'
      IF (scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',U_ref/L_ref,' * $0" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            vort_'//count_char//'.h5:/vort'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
    END IF
    !--- Forces
    IF (write_force_yes) THEN
      !--- X-Force
      WRITE(xmf,fmt1) '        <Attribute Name="Force_X" Center="Node">'
      IF (.NOT. scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',L_ref / (rho_fluid * U_ref**2.),' * $0" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            forceX_'//count_char//'.h5:/forceX'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (.NOT. scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
      !--- Y-Force
      WRITE(xmf,fmt1) '        <Attribute Name="Force_Y" Center="Node">'
      IF (.NOT. scale_output_yes) WRITE(xmf,fmt6) '          <DataItem ItemType="Function" Function="',L_ref / (rho_fluid * U_ref**2.),' * $0" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt3) '          <DataItem Format="HDF" Dimensions="',dim2,dim1,'" NumberType="Float" Precision="8">'
      WRITE(xmf,fmt1) '            forceY_'//count_char//'.h5:/forceY'
      WRITE(xmf,fmt1) '          </DataItem>'
      IF (.NOT. scale_output_yes) WRITE(XMF,fmt1) '          </DataItem>'
      WRITE(xmf,fmt1) '        </Attribute>'
    END IF

    WRITE(xmf,fmt1) '    </Grid>'
    WRITE(xmf,fmt1) '  </Domain>'
    WRITE(xmf,fmt1) '</Xdmf>'

    CLOSE(xmf)

  END SUBROUTINE write_2D_xdmf



  SUBROUTINE write_xdmf_timecollection()

    IMPLICIT NONE

    INTEGER                       ::  xmf=995
    INTEGER                       ::  i
    CHARACTER(LEN=8)              ::  count_char

    CHARACTER(*), PARAMETER       ::  fmt1 = "(a)"

    IF (rank .EQ. 0) THEN
      IF (dimens .EQ. 3) THEN
        IF (write_large) THEN
          OPEN(xmf, FILE='IMPACT_xdmf3d_large.xmf')
          WRITE(xmf,fmt1) '<?xml version="1.0" ?>'
          WRITE(xmf,fmt1) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
          WRITE(xmf,fmt1) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
          WRITE(xmf,fmt1) '  <Domain>'
          WRITE(xmf,fmt1) '    <Grid GridType="Collection" CollectionType="Temporal">'
          DO i= 0, write_count-1
            CALL num_to_string(8, i, count_char)
            WRITE(xmf,fmt1) '      <xi:include href="IMPACT_xdmf3d_large_'//count_char//'.xmf" xpointer="xpointer(/Xdmf/Domain/Grid)" />'
          END DO
          WRITE(xmf,fmt1) '    </Grid>'
          WRITE(xmf,fmt1) '  </Domain>'
          WRITE(xmf,fmt1) '</Xdmf>'
          CLOSE(xmf)
        END IF
        IF (write_med  ) THEN
          OPEN(xmf, FILE='IMPACT_xdmf3d_med.xmf')
          WRITE(xmf,fmt1) '<?xml version="1.0" ?>'
          WRITE(xmf,fmt1) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
          WRITE(xmf,fmt1) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
          WRITE(xmf,fmt1) '  <Domain>'
          WRITE(xmf,fmt1) '    <Grid GridType="Collection" CollectionType="Temporal">'
          DO i= 0, write_count-1
            CALL num_to_string(8, i, count_char)
            WRITE(xmf,fmt1) '      <xi:include href="IMPACT_xdmf3d_med_'//count_char//'.xmf" xpointer="xpointer(/Xdmf/Domain/Grid)" />'
          END DO
          WRITE(xmf,fmt1) '    </Grid>'
          WRITE(xmf,fmt1) '  </Domain>'
          WRITE(xmf,fmt1) '</Xdmf>'
          CLOSE(xmf)
        END IF
        IF (write_small) THEN
          OPEN(xmf, FILE='IMPACT_xdmf3d_small.xmf')
          WRITE(xmf,fmt1) '<?xml version="1.0" ?>'
          WRITE(xmf,fmt1) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
          WRITE(xmf,fmt1) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
          WRITE(xmf,fmt1) '  <Domain>'
          WRITE(xmf,fmt1) '    <Grid GridType="Collection" CollectionType="Temporal">'
          DO i= 0, write_count-1
            CALL num_to_string(8, i, count_char)
            WRITE(xmf,fmt1) '      <xi:include href="IMPACT_xdmf3d_small_'//count_char//'.xmf" xpointer="xpointer(/Xdmf/Domain/Grid)" />'
          END DO
          WRITE(xmf,fmt1) '    </Grid>'
          WRITE(xmf,fmt1) '  </Domain>'
          WRITE(xmf,fmt1) '</Xdmf>'
          CLOSE(xmf)
        END IF
      ELSE
        OPEN(xmf, FILE='IMPACT_xdmf2d.xmf')
        WRITE(xmf,fmt1) '<?xml version="1.0" ?>'
        WRITE(xmf,fmt1) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
        WRITE(xmf,fmt1) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
        WRITE(xmf,fmt1) '  <Domain>'
        WRITE(xmf,fmt1) '    <Grid GridType="Collection" CollectionType="Temporal">'
          DO i= 0, write_count-1
            CALL num_to_string(8, i, count_char)
            WRITE(xmf,fmt1) '      <xi:include href="IMPACT_xdmf2d_'//count_char//'.xmf" xpointer="xpointer(/Xdmf/Domain/Grid)" />'
          END DO
        WRITE(xmf,fmt1) '    </Grid>'
        WRITE(xmf,fmt1) '  </Domain>'
        WRITE(xmf,fmt1) '</Xdmf>'
        CLOSE(xmf)
      END IF
    END IF

  END SUBROUTINE write_xdmf_timecollection


  SUBROUTINE start_hdf5_for_testing

  IMPLICIT NONE

  CALL h5open_f(herror)

  END SUBROUTINE start_hdf5_for_testing




  SUBROUTINE stop_hdf5_for_testing

  IMPLICIT NONE

  CALL h5close_f(herror)

  END SUBROUTINE stop_hdf5_for_testing

END MODULE mod_inout
