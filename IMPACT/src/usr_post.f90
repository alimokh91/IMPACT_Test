!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE postprocess
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE mod_inout
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  
  !--- data post-processing ---
  ! note: - use the following subroutine for post-processing
  !       - the content is just an example/template, feel free to edit and modify it
  !       - currently only subroutines for reading 3D data plus ghost-cell update are available
  !       - available subroutines for 2D can only read (without ghost-cell update)
  !
  INTEGER                ::  m
  CHARACTER(LEN=5)       ::  count_char
  CHARACTER(LEN=1)       ::  conc_number
  
  
                         OPEN(21,FILE='post_TKE_mod_restart'       //restart_char//'.txt',STATUS='UNKNOWN')
                         OPEN(22,FILE='post_TKE_spectrum_restart'  //restart_char//'.txt',STATUS='UNKNOWN')
                         OPEN(23,FILE='post_turb_stat_restart'     //restart_char//'.txt',STATUS='UNKNOWN')
  
  
  DO timestep = 1, 50
     CALL num_to_string(5,timestep,count_char)
     CALL read2_hdf('velX_med_'//count_char,'velX' ,S1p,S2p,S3p,N1p,N2p,N3p,0,work1)
     CALL read2_hdf('velY_med_'//count_char,'velY' ,S1p,S2p,S3p,N1p,N2p,N3p,0,work2)
     CALL read2_hdf('velZ_med_'//count_char,'velZ' ,S1p,S2p,S3p,N1p,N2p,N3p,0,work3)

     CALL compute_stats
  END DO
  
  CALL close_stats
  
  
  END SUBROUTINE postprocess
