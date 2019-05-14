!*************************************************************************************************************
!* IB_IMPACT  												     *
!* by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)      			     *
!* January 2015                                                                                              *
!*************************************************************************************************************


MODULE mod_ibm

  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_lib !num_to_string
  !USE mpi

  USE omp_lib
  

  PRIVATE

  PUBLIC init_ib_mpi
  PUBLIC init_ib
  PUBLIC spread_force_dens_to_vel_grid
  PUBLIC interpolate_vel_to_ib
  PUBLIC find_vicinity_coordinates
  PUBLIC update_boundary
  PUBLIC calculate_displacements
  PUBLIC discrete_delta
  PUBLIC compute_strain_stress
  PUBLIC compute_force
  PUBLIC calculate_node_volumes
  PUBLIC fe_setup_triangular_3node
  PUBLIC write_hdf_ib
  PUBLIC setup_mass_matrix
  PUBLIC setup_damping_matrix
  PUBLIC setup_stiffness_matrix
  PUBLIC gq_int
  PUBLIC qp_lookup
  PUBLIC integrand

  INCLUDE 'mpif.h'

  CONTAINS

  SUBROUTINE init_ib_mpi

  IMPLICIT NONE

  INTEGER                :: world_size, myrank, color

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, merror)

  IF (myrank .EQ. (world_size-1)) THEN
    color = 2
  ELSE
    color = 1
  END IF

  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, color, 1, COMM_LOCAL, merror)
  CALL MPI_INTERCOMM_CREATE(COMM_LOCAL,0,MPI_COMM_WORLD,0,111,COMM_INTER,merror)

  END SUBROUTINE init_ib_mpi
  
  !> subroutine to initialize the IB, setup parameters, allocate fields and to make sure only one process is assigend with task
  SUBROUTINE init_ib
  
  IMPLICIT NONE

  INTEGER                :: i, j, ios
  REAL                   :: pi
  REAL                   :: delta_b
  CHARACTER(len=80)      :: text, dummy
  INTEGER, ALLOCATABLE   :: dof_red(:)
  INTEGER                :: dof_red_len
  REAL, ALLOCATABLE      :: m_mass_red(:,:), k_stiff_red(:,:)
  REAL, ALLOCATABLE      :: alpha_r(:), alpha_i(:), beta_r(:), work_ev(:)
  LOGICAL, ALLOCATABLE   :: mask(:,:)
  REAL                   :: DUMMY_EV(1,1)
  INTEGER                :: ev_info
  !--- Oscillator Model aorta ---
  INTEGER                :: ii, jj
  INTEGER                :: n_points_lev, n_points_z, n_points_total
  REAL                   :: delta_z, l_half, phi_aorta, theta_aorta

  pi = 2.*ABS(ACOS(0.))

  IF (fem_yes .AND. dimens .EQ. 2) THEN
  
  OPEN(20,FILE=TRIM('boundary.txt'),FORM='formatted',ACTION='read',STATUS='old',IOSTAT=ios)

  IF (ios /=0) THEN
    IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open boundary.txt file!'
    CALL MPI_FINALIZE(merror)
    STOP
  END IF

  !=== read in number of nodes and elements and BCs ===========================================================
  DO 
    READ(20,FMT='(a)',IOSTAT=ios) text
   
    !write(*,*) text

    IF (ios /=0) THEN
      IF (rank == 0) WRITE(*,*) 'ERROR! Cannot read boundary.txt file!'
      CALL MPI_FINALIZE(merror)
      STOP
    END IF

  !--- boundary nodes -----------------------------------------------------------------------------------------
    IF (INDEX(text,'nodestotal') == 1) READ(UNIT=text,FMT=*) dummy, M_bound
  !--- elements -----------------------------------------------------------------------------------------------
    IF (INDEX(text,'elemstotal') == 1) READ(UNIT=text,FMT=*) dummy, M_elems
  !--- bcs ----------------------------------------------------------------------------------------------------
    IF (INDEX(text,'bcstotal') == 1) READ(UNIT=text,FMT=*) dummy, M_ebcs

    IF (INDEX(text,'FILEEND') == 1) EXIT
  END DO

  ELSE IF (.NOT. fem_yes .AND. dimens .EQ. 3) THEN
    ! half the length of the secant ~ half the length of the optimal delta_l
    l_half     = reach/4.
    ! delta_l along the aorta
    delta_z    = reach/2.
    n_points_z = NINT((aorta_end - aorta_start)/delta_z) + 1
    
    ! angle of radius to secant
    theta_aorta = ACOS(l_half/aorta_radius)
    ! angle of secant
    phi_aorta   = pi - 2*theta_aorta

    ! number of points on a level
    n_points_lev = NINT(2*pi/phi_aorta)
    ! correct phi so the points are equidistant along the periphery
    phi_aorta   = 2*pi/n_points_lev
    n_points_total = n_points_lev * n_points_z

    M_bound = n_points_total
    M_ebcs  = 3*n_points_lev

  END IF ! fem_yes

  !*** deprecated *********************************************************************************************
  !=== define parameters ======================================================================================
  !--- number of boundary nodes --- (must be even)
  !M_bound = 192

  !--- number of finite elements ---
  !M_elems = 2*(2*M_bound/3 - 2
  !************************************************************************************************************

  !=== allocate fields ========================================================================================
  !--- displacements/Lagrangian values ---
  ALLOCATE(yb(1:M_bound, 1:dimens))
  ALLOCATE(xb(1:M_bound, 1:dimens))
  ALLOCATE(ub(1:M_bound, 1:dimens))
  ALLOCATE(db(1:M_bound, 1:dimens))
  ALLOCATE(container(1:M_bound, 1:dimens))

  !--- Force fields ---
  ALLOCATE(fb(1:M_bound, 1:dimens))
  ! bbecsek 110216: Moved to 'alloc.f90'
  !ALLOCATE(fd(b1L:(N1+b1U), b2L:(N2+b2U), b3L:(N3+b3U),1:3)) 
  
  !--- Discretization area/volume elements ---
  ALLOCATE(node_vol(1:M_bound))

  !--- Zero displacement bcs ---
  ALLOCATE(ebcs(1:M_ebcs,1:(dimens+1)))


  IF (fem_yes .AND. dimens .EQ. 2) THEN

  !--- FE variables ---
  ALLOCATE(elems(1:M_elems, 1:3))
  ALLOCATE(strain(1:M_elems, 1:3))
  ALLOCATE(stress(1:M_elems, 1:3))
  ALLOCATE(k_stiff(1:2*M_bound, 1:2*M_bound))
  ALLOCATE(m_mass (1:2*M_bound, 1:2*M_bound))
  ALLOCATE(c_damp (1:2*M_bound, 1:2*M_bound))

  !=== MPI ====================================================================================================
  !--- gather total number of involved processes ---
  CALL MPI_COMM_SIZE(COMM_CART,total_ranks,merror)

  END IF ! fem_yes

  !=== initialize fields ======================================================================================
  !--- displacements/Lagrangian coordinates ---
  yb          = 0.
  xb          = 0.
  ub          = 0.
  db          = 0.
  container   = 0.

  !--- Force fields ---
  fb          = 0. 
  fd          = 0.

  !--- Discretization area/volume elements ---
  node_vol    = 0.
  
  !--- Zero displacement bcs ---
  ebcs        = 0.

  IF (fem_yes) THEN
  !--- FE variables ---
  elems       = 0
  strain      = 0.
  stress      = 0.
  k_stiff     = 0.
  ebcs        = 0.
  
  !=== initialize boundary ====================================================================================
  !*** deprecated *********************************************************************************************
  !delta_b = L2/(2.*(REAL(M_bound)/3. - 1.))


  !delta_b1 = delta_b
  !delta_b2 = delta_b
  !delta_b3 = delta_b

  !DO i = 1, M_bound
  !  IF (i .LE. M_bound/3) THEN
  !    yb(i,2) = REAL(i-1)*(delta_b)
  !    yb(i,1) = L1 - 2.5
  !  ELSE IF (i .GT. M_bound/3 .AND. i .LE. 2*M_bound/3) THEN
  !    yb(i,2) = REAL(2*M_bound/3-i)*(delta_b)
  !    yb(i,1) = L1 - 2.5 + (delta_b)
  !  ELSE IF (i .GT. 2*M_bound/3) THEN
  !    yb(i,2) = REAL(i-1-2*M_bound/3)*(delta_b)
  !    yb(i,1) = L1 - 2.5 + 2.*(delta_b)
  !  END IF
  !END DO
  !************************************************************************************************************
  
  REWIND(20,IOSTAT=ios)

  !=== read nodes from file ===================================================================================
  DO
    READ(20,FMT=*,IOSTAT=ios) text
    IF (text == 'NODESBEGIN') EXIT
  END DO

  DO
    READ(20,FMT='(a)',IOSTAT=ios) text
    IF (INDEX(text,'NODESEND') == 1) EXIT
    READ(UNIT=text,FMT=*) i, yb(i,1:2)
    !IF(rank==0) WRITE(*,*) 'Reading node ', i,' of ', M_bound, ' coordinates...'
  END DO

  !yb(:,1) = yb(:,1) + L1-2.5 !< @note: rather add this directly to the file boundary.txt
  
  !=== read element nodes from file ===========================================================================
  DO
    READ(20,FMT=*,IOSTAT=ios) text
    IF (text == 'ELEMSBEGIN') EXIT
  END DO

  DO 
    READ(20,FMT='(a)',IOSTAT=ios) text
    IF (INDEX(text,'ELEMSEND') == 1) EXIT
    READ(UNIT=text,FMT=*) i, elems(i,1:3)
    !IF(rank==0) WRITE(*,*) 'Reading element ', i ,' of ', M_elems, ' nodes ...'
  END DO

  !=== read essential BCs from file ===========================================================================
  DO
    READ(20,FMT=*,IOSTAT=ios) text
    IF (text == 'BCSBEGIN') EXIT
  END DO

  DO i = 1, M_ebcs
    READ(20,FMT='(a)',IOSTAT=ios) text
    IF (INDEX(text,'BCSEND') == 1) EXIT
    READ(UNIT=text,FMT=*) ebcs(i,1:3)
    !IF(rank==0) WRITE(*,*) 'Reading displacement BCs for node ', NINT(ebcs(i,1))
  END DO

  CLOSE(20)

  ELSE IF (.NOT. fem_yes .AND. dimens .EQ. 3) THEN
    DO i = 1, M_bound
      ii = MOD(i, n_points_lev)
      jj = (i-1)/n_points_lev ! floating point arithmatics
      yb(i,1) = jj*delta_z + aorta_start
      yb(i,2) = aorta_radius*COS(ii*phi_aorta) + 1.
      yb(i,3) = aorta_radius*SIN(ii*phi_aorta) + 1.
      node_vol(i) = delta_z**3.
    END DO
    DO i = 1, 3*n_points_lev
      !ebcs(i               ,1) = i
      !ebcs(i+  n_points_lev,1) = i + n_points_lev
      !ebcs(i+2*n_points_lev,1) = n_points_total - i + 1
      ebcs(i,1) = -10
    END DO
  END IF ! fem_yes

  !--- save reference configuration ---
  xb = yb
  db = yb - xb

  !*** deprecated *********************************************************************************************
  !=== define finite elements (2D triangular 3-node) ==========================================================
  !--- Make sure element nodes are always defined in counter-clockwise order! ---
  !DO i = 1, M_elems
  !  IF (i .LE. (2*M_bound/3 -2)) THEN
  !    IF (MOD(i,2) .EQ. 0) THEN
  !      elems(i,1) = i/2
  !      elems(i,2) = 2*M_bound/3 - i/2
  !      elems(i,3) = i/2 + 1
  !    ELSE
  !      elems(i,1) = (i+1)/2 
  !      elems(i,2) = 2*M_bound/3 - (i-1)/2
  !      elems(i,3) = 2*M_bound/3 - (i+1)/2
  !    END IF
  !  ELSE IF (i .GT. (2*M_bound/3-2)) THEN
  !    IF (MOD(i,2) .EQ. 0) THEN
  !      elems(i,1) = 2*M_bound/3 - (i - 2*M_bound/3)/2
  !      elems(i,2) = 2*M_bound/3 + (i - 2*M_bound/3 + 4)/2
  !      elems(i,3) = 2*M_bound/3 - (i - 2*M_bound/3 + 2)/2
  !    ELSE
  !      elems(i,1) = 2*M_bound/3 - (i - 2*M_bound/3 + 1)/2
  !      elems(i,2) = 2*M_bound/3 + (i - 2*M_bound/3 + 3)/2
  !      elems(i,3) = 2*M_bound/3 + (i - 2*M_bound/3 + 5)/2
  !    END IF
  !  END IF
  !  !if(rank==0)write(*,*) 'elem',i,'a:',elems(i,1),'b:',elems(i,2),'c:',elems(i,3)
  !END DO
    !**********************************************************************************************************

  IF (fem_yes .AND. dimens .EQ. 2) THEN

  !=== calculate nodes vol's in initial state =================================================================
  DO i = 1, M_elems
    CALL calculate_node_volumes(i)
  END DO

  !=== define strain-stress relationship (Hooke's Law) ========================================================
  C(1,1) = 1.
  C(1,2) = nu_poiss
  C(1,3) = 0.
  C(2,1) = nu_poiss
  C(2,2) = 1.
  C(2,3) = 0.
  C(3,1) = 0.
  C(3,2) = 0.
  C(3,3) = (1.-nu_poiss)/2.

  C = (E_mod/(1.-nu_poiss**2.))*C
 
  !=== setup global element matrices ==========================================================================
  CALL fe_setup_triangular_3node

  !=== do vibration free analysis to determine eigenfrequencies ===============================================
  dof_red_len = SIZE(ebcs,1)
  ALLOCATE(dof_red(1:(2*dof_red_len)))
  ALLOCATE(m_mass_red(1:(2*(M_bound-dof_red_len)),1:(2*(M_bound-dof_red_len))))
  ALLOCATE(k_stiff_red(1:(2*(M_bound-dof_red_len)),1:(2*(M_bound-dof_red_len))))
  ALLOCATE(mask(1:(2*M_bound),1:(2*M_bound)))
  mask = .True.
  ALLOCATE(alpha_r(1:2*(M_bound-dof_red_len)))
  ALLOCATE(alpha_i(1:2*(M_bound-dof_red_len)))
  ALLOCATE(beta_r(1:2*(M_bound-dof_red_len)))
  ALLOCATE(work_ev(1:16*(M_bound-dof_red_len)))

  !--- determine DOFs with essential BCs ---
  DO i = 1, dof_red_len
    dof_red(2*i-1) = NINT(ebcs(i,1)*2 - 1)
    dof_red(2*i)   = NINT(ebcs(i,1)*2)
  END DO

  !--- create mask for deleting matrix row/columns ---
  DO i = 1, 2*M_bound
    IF ( ANY( i == dof_red ) ) mask(i,:) = .False.
    DO j = 1, 2*M_bound
      IF (ANY( j == dof_red ) ) mask(:,j) = .False.
    END DO
  END DO

  !--- delete matrix rows/columns ---
  m_mass_red  = RESHAPE(PACK(m_mass, mask), (/2*(M_bound-dof_red_len), 2*(M_bound-dof_red_len) /) )
  k_stiff_red = RESHAPE(PACK(k_stiff,mask), (/2*(M_bound-dof_red_len), 2*(M_bound-dof_red_len) /) )
  m_mass_red  = TRANSPOSE(m_mass_red)
  k_stiff_red = TRANSPOSE(k_stiff_red)

  !--- solve EV-Problem to obtain eigenfrequencies ---
  IF (rank==0) WRITE(*,*) 'Solving structural eigenvalue problem...'
  CALL DGGEV('N','N',2*(M_bound-dof_red_len),k_stiff_red,2*(M_bound-dof_red_len),m_mass_red,2*(M_bound-dof_red_len), & 
         alpha_r, alpha_i, beta_r, DUMMY_EV, 1, DUMMY_EV, 1, work_ev, 16*(M_bound-dof_red_len), ev_info)

  !--- check outcome ---
  IF (ev_info .GT. 0) THEN
    WRITE(*,*) 'WARNING: Failure in DGGEV. INFO = ', ev_info
  ELSE
    max_ev_freq = 0.
    DO i = 1, 2*(M_bound - dof_red_len)
      IF (alpha_i(i) .NE. 0.0D0) THEN
        WRITE(*,*) 'There is at least one complex EV.'
        IF (beta_r(i) .NE. 0.0D0) max_ev_freq = MAX(max_ev_freq, SQRT(ABS(alpha_r(i))/ABS(beta_r(i))))
      ELSE
        IF (beta_r(i) .NE. 0.0D0) max_ev_freq = MAX(max_ev_freq, SQRT(    alpha_r(i) /    beta_r(i) ))
      END IF
    END DO
  END IF

  IF (rank==0) WRITE(*,*)'   ...fastest mode = ', max_ev_freq

  !--- free memory ---
  DEALLOCATE(dof_red)
  DEALLOCATE(m_mass_red)
  DEALLOCATE(k_stiff_red)
  DEALLOCATE(mask)
  DEALLOCATE(alpha_r)
  DEALLOCATE(alpha_i)
  DEALLOCATE(beta_r)
  DEALLOCATE(work_ev)

  END IF ! fem_yes

  END SUBROUTINE init_ib


  

  !> subroutine to spread the individual force components from the Lagrangian grid onto the individual velocity grids
  !! @todo: make it compatible with non-equidistant meshes
  SUBROUTINE spread_force_dens_to_vel_grid
    
    IMPLICIT NONE 
    
    INTEGER                             :: i, j, k, l
    INTEGER, DIMENSION(2)               :: ind_1p, ind_2p, ind_3p, ind_1u, ind_2v, ind_3w
    REAL   , DIMENSION(:), ALLOCATABLE  :: dist_1p, dist_2p, dist_3p, dist_1u, dist_2v, dist_3w
    INTEGER, DIMENSION(shape(req_bcast)):: status(MPI_STATUS_SIZE)

    
    
   
    !--- first thing's first: empty the force density container before adding anything ---
    fd = 0.

    CALL MPI_BARRIER(COMM_CART,merror)
    !--- broadcast boundary force and coordinates to all processe under ---
    !--- communicator COMM_CART ---
    CALL MPI_IBCAST(fb,SIZE(fb),MPI_REAL8,0,COMM_CART,req_bcast(1),merror)
    CALL MPI_IBCAST(yb,SIZE(yb),MPI_REAL8,0,COMM_CART,req_bcast(2),merror)
    CALL MPI_IBCAST(node_vol,SIZE(node_vol),MPI_REAL8,0,COMM_CART,req_bcast(3),merror)
    CALL MPI_WAITALL(3,req_bcast,status,merror)
    !--- loop over all Lagrangian points ---
    SPREAD_LOOP: DO i = 1, M_bound
      IF ( ANY( i == NINT(ebcs(:,1) ) ) ) cycle SPREAD_LOOP
      !--- first, find all Eulerian coordinates that are within 'reach' from the Lagrangian point ---
      IF (dimens .EQ. 3) THEN
        CALL find_vicinity_coordinates(i, ind_1p, ind_1u, ind_2p, ind_2v, ind_3p, ind_3w)
      ELSE
        CALL find_vicinity_coordinates(i, ind_1p, ind_1u, ind_2p, ind_2v)
      END IF
      !--- allocate fields for distances ---
      ALLOCATE( dist_1p(1:(ind_1p(2)-ind_1p(1) + 1)) )
      ALLOCATE( dist_2p(1:(ind_2p(2)-ind_2p(1) + 1)) )
      ALLOCATE( dist_1u(1:(ind_1u(2)-ind_1u(1) + 1)) )
      ALLOCATE( dist_2v(1:(ind_2v(2)-ind_2v(1) + 1)) )
      IF (dimens == 3) THEN
      ALLOCATE( dist_3p(1:(ind_3p(2)-ind_3p(1) + 1)) )
      ALLOCATE( dist_3w(1:(ind_3w(2)-ind_3w(1) + 1)) )
      END IF
      !--- second, compute the distance of these coordinates to the Lagrangian location ---
      IF (dimens .EQ. 3) THEN
        CALL compute_distances_to_vicinity_coordinates(i, ind_1p , ind_1u , ind_2p , ind_2v , ind_3p , ind_3w , &
                                                          dist_1p, dist_1u, dist_2p, dist_2v, dist_3p, dist_3w)
      ELSE
        CALL compute_distances_to_vicinity_coordinates(i, ind_1p , ind_1u , ind_2p , ind_2v , (/0, 0/) , (/0, 0/) , &
                                                          dist_1p, dist_1u, dist_2p, dist_2v)
      END IF

      !*** debugging
      !write(*,*) 'in process rank=',rank,' ind_1p=',ind_1p
      !write(*,*) 'in process rank=',rank,' ind_2p=',ind_2p
      !write(*,*) 'in process rank=',rank,' ind_3p=',ind_3p
      !write(*,*) 'in process rank=',rank,' ind_1u=',ind_1u
      !write(*,*) 'in process rank=',rank,' ind_2v=',ind_2v
      !write(*,*) 'in process rank=',rank,' ind_3w=',ind_3w
      !write(*,*) 'in process rank=',rank,' dist_1p=',dist_1p
      !write(*,*) 'in process rank=',rank,' dist_2p=',dist_2p
      !write(*,*) 'in process rank=',rank,' dist_3p=',dist_3p
      !write(*,*) 'in process rank=',rank,' dist_1u=',dist_1u
      !write(*,*) 'in process rank=',rank,' dist_2v=',dist_2v
      !write(*,*) 'in process rank=',rank,' dist_3w=',dist_3w
      !if(rank==0)write(*,*)  'in process rank=',rank,' fb =',fb

      !--- when done, loop over the indices and spread the force over the Eulerian grid, vel-comp-wise ---
      IF (dimens == 3) THEN
      !=== vel comp 1 =========================================================================================
        DO l = ind_3p(1), ind_3p(2)
          DO k = ind_2p(1), ind_2p(2)
            DO j = ind_1u(1), ind_1u(2)
              fd(j,k,l,1) = fd(j,k,l,1) + fb(i,1) * discrete_delta(dist_1u( (j - ind_1u(1) + 1) ))&
                                                          * discrete_delta(dist_2p( (k - ind_2p(1) + 1) ))&
                                                          * discrete_delta(dist_3p( (l - ind_3p(1) + 1) ))&
                                                          * node_vol(i) ! this should be superfluous, since we need to divide by some 'volume' to get a force density
            END DO
          END DO
        END DO
      !=== vel comp 2 =========================================================================================
        DO l = ind_3p(1), ind_3p(2)
          DO k = ind_2v(1), ind_2v(2)
            DO j = ind_1p(1), ind_1p(2)
              fd(j,k,l,2) = fd(j,k,l,2) + fb(i,2) * discrete_delta(dist_1p( (j - ind_1p(1) + 1) ))&
                                                          * discrete_delta(dist_2v( (k - ind_2v(1) + 1) ))&
                                                          * discrete_delta(dist_3p( (l - ind_3p(1) + 1) ))&
                                                          * node_vol(i) ! see above
            END DO
          END DO
        END DO
      !=== vel comp 3 =========================================================================================
        DO l = ind_3w(1), ind_3w(2)
          DO k = ind_2p(1), ind_2p(2)
            DO j = ind_1p(1), ind_1p(2)
              fd(j,k,l,3) = fd(j,k,l,3) + fb(i,3) * discrete_delta(dist_1p( (j - ind_1p(1) + 1) ))&
                                                          * discrete_delta(dist_2p( (k - ind_2p(1) + 1) ))&
                                                          * discrete_delta(dist_3w( (l - ind_3w(1) + 1) ))&
                                                          * node_vol(i) ! see above
            END DO
          END DO
        END DO
        
      ELSE
      !=== vel comp 1 =========================================================================================
        l = 1
          DO k = ind_2p(1), ind_2p(2)
            DO j = ind_1u(1), ind_1u(2)
              fd(j,k,l,1) = fd(j,k,l,1) + fb(i,1) * discrete_delta(dist_1u( (j - ind_1u(1) + 1) ))&
                                                          * discrete_delta(dist_2p( (k - ind_2p(1) + 1) ))&
                                                          * node_vol(i) ! see above
            END DO
          END DO
        
      !=== vel comp 2 =========================================================================================
        l = 1
          DO k = ind_2v(1), ind_2v(2)
            DO j = ind_1p(1), ind_1p(2)
              fd(j,k,l,2) = fd(j,k,l,2) + fb(i,2) * discrete_delta(dist_1p( (j - ind_1p(1) + 1) ))&
                                                          * discrete_delta(dist_2v( (k - ind_2v(1) + 1) ))&
                                                          * node_vol(i) ! see above
            END DO
          END DO
        
      END IF
      !--- deallocate fields for distances ---
      DEALLOCATE( dist_1p )
      DEALLOCATE( dist_2p )
      DEALLOCATE( dist_1u )
      DEALLOCATE( dist_2v )
      IF (dimens == 3) THEN
      DEALLOCATE( dist_3p )
      DEALLOCATE( dist_3w )
      END IF


    END DO SPREAD_LOOP

   !write(*,*)'in process',rank,'sum_fd=',sum(fd)


    !--- exchange the force field throughout the processes ---
    CALL exchange_part(1,1,1,fd)
    CALL exchange_part(2,1,1,fd)
    CALL exchange_part(3,1,1,fd)

    CALL exchange_part(1,2,1,fd)
    CALL exchange_part(2,2,1,fd)
    CALL exchange_part(3,2,1,fd)

    IF (dimens .EQ. 3) THEN
      CALL exchange_part(1,3,1,fd)
      CALL exchange_part(2,3,1,fd)
      CALL exchange_part(3,3,1,fd)
    END IF

   !write(*,*)'in process',rank,'sum_fd=',sum(fd)

  END SUBROUTINE spread_force_dens_to_vel_grid




  !> subroutine to interpolate the individual velocity components onto the Lagrangian grid
  SUBROUTINE interpolate_vel_to_ib

    IMPLICIT NONE
    
    INTEGER                            :: i, j, k, l
    INTEGER, DIMENSION(2)              :: ind_1p, ind_2p, ind_3p, ind_1u, ind_2v, ind_3w
    REAL   , DIMENSION(:), ALLOCATABLE :: dist_1p, dist_2p, dist_3p, dist_1u, dist_2v, dist_3w
    REAL                               :: u1_reduced, u2_reduced, u3_reduced, u1_forsend, u2_forsend, u3_forsend
    INTEGER, DIMENSION(shape(req_red)) :: status(MPI_STATUS_SIZE)

    !--- store old entries and delete any previous ones of ub ---
    container = 0.
    container = ub
    ub        = 0.

    CALL MPI_BARRIER(COMM_CART,merror)
    !--- broadcast boundary coordinate across all processes ---
    !CALL MPI_IBCAST(yb,SIZE(yb),MPI_REAL8,0,COMM_CART,req_bcast,merror)
    !CALL MPI_WAIT(req_bcast,status,merror)
    !if(rank==0)write(*,*)  'in process rank=',rank,' ubb =',yb

    !--- loop over all Lagrangian points ---
    INTERP_LOOP: DO i = 1, M_bound
      IF ( ANY( i == NINT(ebcs(:,1) ) ) ) cycle INTERP_LOOP
      !--- first, find all Eulerian coordinates that are within 'reach' from the Lagrangian point ---
      IF (dimens .EQ. 3) THEN
        CALL find_vicinity_coordinates(i, ind_1p, ind_1u, ind_2p, ind_2v, ind_3p, ind_3w)
      ELSE
        CALL find_vicinity_coordinates(i, ind_1p, ind_1u, ind_2p, ind_2v)
      END IF
      !--- allocate fields for distances ---
      ALLOCATE( dist_1p(1:(ind_1p(2)-ind_1p(1) + 1)) )
      ALLOCATE( dist_2p(1:(ind_2p(2)-ind_2p(1) + 1)) )
      ALLOCATE( dist_1u(1:(ind_1u(2)-ind_1u(1) + 1)) )
      ALLOCATE( dist_2v(1:(ind_2v(2)-ind_2v(1) + 1)) )
      IF (dimens == 3) THEN
      ALLOCATE( dist_3p(1:(ind_3p(2)-ind_3p(1) + 1)) )
      ALLOCATE( dist_3w(1:(ind_3w(2)-ind_3w(1) + 1)) )
      END IF
      !--- second, compute the distances of these coordinates to the Lagrangian location ---
      IF (dimens .EQ. 3) THEN
        CALL compute_distances_to_vicinity_coordinates(i, ind_1p , ind_1u , ind_2p , ind_2v , ind_3p , ind_3w , &
                                                          dist_1p, dist_1u, dist_2p, dist_2v, dist_3p, dist_3w)
      ELSE
        CALL compute_distances_to_vicinity_coordinates(i, ind_1p , ind_1u , ind_2p , ind_2v , (/0, 0/) , (/0, 0/) , &
                                                          dist_1p, dist_1u, dist_2p, dist_2v)
      END IF
      !--- when done, loop over the indices and sum them up to equal the Lagrangian velocity, vel-comp-wise ---
      IF (dimens == 3) THEN
      !=== vel comp 1 =========================================================================================
        DO l = ind_3p(1), ind_3p(2)
          DO k = ind_2p(1), ind_2p(2)
            DO j = ind_1u(1), ind_1u(2)
              ub(i, 1) = ub(i, 1) + discrete_delta(dist_1u( (j - ind_1u(1) + 1) ))&
                                  * discrete_delta(dist_2p( (k - ind_2p(1) + 1) ))&
                                  * discrete_delta(dist_3p( (l - ind_3p(1) + 1) ))&
                                  * vel(j,k,l,1)*dx1u(j)*dx2p(k)*dx3p(l)
            END DO
          END DO
        END DO
      !=== vel comp 2 =========================================================================================
        DO l = ind_3p(1), ind_3p(2)
          DO k = ind_2v(1), ind_2v(2)
            DO j = ind_1p(1), ind_1p(2)
              ub(i, 2) = ub(i, 2) + discrete_delta(dist_1p( (j - ind_1p(1) + 1) ))&
                                  * discrete_delta(dist_2v( (k - ind_2v(1) + 1) ))&
                                  * discrete_delta(dist_3p( (l - ind_3p(1) + 1) ))&
                                  * vel(j,k,l,2)*dx1p(j)*dx2v(k)*dx3p(l)
            END DO
          END DO
        END DO
      !=== vel comp 3 =========================================================================================
        DO l = ind_3w(1), ind_3w(2)
          DO k = ind_2p(1), ind_2p(2)
            DO j = ind_1p(1), ind_1p(2)
              ub(i, 3) = ub(i, 3) + discrete_delta(dist_1p( (j - ind_1p(1) + 1) ))&
                                  * discrete_delta(dist_2p( (k - ind_2p(1) + 1) ))&
                                  * discrete_delta(dist_3w( (l - ind_3w(1) + 1) ))&
                                  * vel(j,k,l,3)*dx1p(j)*dx2p(k)*dx3w(l)
            END DO
          END DO
        END DO
      ELSE
      !=== vel comp 1 =========================================================================================
        l = 1
          DO k = ind_2p(1), ind_2p(2)
            DO j = ind_1u(1), ind_1u(2)
              ub(i, 1) = ub(i, 1) + discrete_delta(dist_1u( (j - ind_1u(1) + 1) ))&
                                  * discrete_delta(dist_2p( (k - ind_2p(1) + 1) ))&
                                  * vel(j,k,l,1)*dx1u(j)*dx2p(k)
            END DO
          END DO
         
      !=== vel comp 2 =========================================================================================
        l = 1
          DO k = ind_2v(1), ind_2v(2)
            DO j = ind_1p(1), ind_1p(2)
              ub(i, 2) = ub(i, 2) + discrete_delta(dist_1p( (j - ind_1p(1) + 1) ))&
                                  * discrete_delta(dist_2v( (k - ind_2v(1) + 1) ))&
                                  * vel(j,k,l,2)*dx1p(j)*dx2v(k)
            END DO
          END DO
        
      END IF
      !--- deallocate fields for distances ---
      DEALLOCATE( dist_1p )
      DEALLOCATE( dist_2p )
      DEALLOCATE( dist_1u )
      DEALLOCATE( dist_2v )
      IF (dimens == 3) THEN
      DEALLOCATE( dist_3p )
      DEALLOCATE( dist_3w )
      END IF



    END DO INTERP_LOOP

    !write(*,*) 'process',rank,'sum_ub=',sum(ub)

    CALL MPI_BARRIER(COMM_CART,merror)
    DO i=1,M_bound
      !--- reduce the boundary velocity to root process ---
                         u1_forsend = ub(i,1)
                         u2_forsend = ub(i,2)
      IF (dimens .EQ. 3) u3_forsend = ub(i,3)
                         CALL MPI_IREDUCE(u1_forsend,u1_reduced,1,MPI_REAL8,MPI_SUM,0,COMM_CART,req_red(1),merror)
                         CALL MPI_IREDUCE(u2_forsend,u2_reduced,1,MPI_REAL8,MPI_SUM,0,COMM_CART,req_red(2),merror)
      IF (dimens .EQ. 3) CALL MPI_IREDUCE(u3_forsend,u3_reduced,1,MPI_REAL8,MPI_SUM,0,COMM_CART,req_red(3),merror)
      
      IF (dimens .NE. 3) THEN
        CALL MPI_WAITALL(2,req_red(1:2),status,merror)
      ELSE
        CALL MPI_WAITALL(3,req_red(1:3),status,merror)
      ENDIF

      IF(rank == 0) THEN
                         ub(i,1) = u1_reduced
                         ub(i,2) = u2_reduced
      IF (dimens .EQ. 3) ub(i,3) = u3_reduced
      END IF
    END DO
    !if(rank==0)write(*,*)  'in process rank=',rank,' sum_ubold =', sum(container)
    !if(rank==0)write(*,*)  'in process rank=',rank,' sum_ub =', sum(ub)

 


  END SUBROUTINE interpolate_vel_to_ib




  !> subroutine to find all Eulerian coordinates' indices that are within the reach 'reach' of the kernel functions
  SUBROUTINE find_vicinity_coordinates(Lag_ind, Eul_ind_1p, Eul_ind_1u, Eul_ind_2p, Eul_ind_2v, Eul_ind_3p, Eul_ind_3w)

  IMPLICIT NONE

  INTEGER, INTENT(IN)                          :: Lag_ind
  INTEGER, INTENT(OUT), DIMENSION(2)           :: Eul_ind_1p, Eul_ind_2p, Eul_ind_1u, Eul_ind_2v
  INTEGER, INTENT(OUT), DIMENSION(2), OPTIONAL :: Eul_ind_3p, Eul_ind_3w

  INTEGER                            :: i, j 
  
 
  !--- compare the Lagrangian coordinate's first entry against the x1p grid (for vel's 2 and 3) ---
  FIND_1p: DO i = S1p, N1p
    !--- add the first and last Eulerian coordinate's index that is within 'reach' to the list ---
    IF ( abs(yb(Lag_ind,1) - x1p(i)) .LE. reach ) THEN
      Eul_ind_1p(1) = i
      DO j = Eul_ind_1p(1), N1p
        IF ( (x1p(j) - yb(Lag_ind,1)) .GT. reach ) THEN
          Eul_ind_1p(2) = (j-1)
          EXIT FIND_1p
        ELSE IF ( j == N1p ) THEN
          Eul_ind_1p(2) = j
          EXIT FIND_1p
        END IF
      END DO
    ELSE
      Eul_ind_1p = -1
    END IF
  END DO FIND_1p

  !--- compare the Lagrangian coordinate's second entry against the x2p grid (for vel's 1 and 3) ---
  FIND_2p: DO i = S2p, N2p
    !--- add the first and last Eulerian coordinate's index that is within 'reach' to the list ---
    IF ( abs(yb(Lag_ind,2) - x2p(i)) .LE. reach ) THEN
      Eul_ind_2p(1) = i
      DO j = Eul_ind_2p(1), N2p
        IF ( (x2p(j) - yb(Lag_ind,2)) .GT. reach ) THEN
          Eul_ind_2p(2) = (j-1)
          EXIT FIND_2p
        ELSE IF ( j == N2p ) THEN
          Eul_ind_2p(2) = j
          EXIT FIND_2p
        END IF
      END DO
    ELSE
      Eul_ind_2p = -1
    END IF
  END DO FIND_2p

  IF (dimens == 3) THEN
  !--- compare the Lagrangian coordinate's third entry against the x3p grid (for vel's 1 and 2) ---
  FIND_3p: DO i = S3p, N3p
    !--- add the first and last Eulerian coordinate's index that is within 'reach' to the list ---
    IF ( abs(yb(Lag_ind,3) - x3p(i)) .LE. reach ) THEN
      Eul_ind_3p(1) = i
      DO j = Eul_ind_3p(1), N3p
        IF ( (x3p(j) - yb(Lag_ind,3)) .GT. reach ) THEN
          Eul_ind_3p(2) = (j-1)
          EXIT FIND_3p
        ELSE IF ( j == N3p ) THEN 
          Eul_ind_3p(2) = j
          EXIT FIND_3p
        END IF
      END DO
    ELSE
      Eul_ind_3p = -1
    END IF
  END DO FIND_3p
  END IF

  !--- compare the Lagrangian coordinate's first entry against the x1u grid (for vel 1) ---
  FIND_1u: DO i = S11, N11 !S11B, N11B?
    !--- add the first and last Eulerian coordinate's index that is within 'reach' to the list ---
    IF ( abs(yb(Lag_ind,1) - x1u(i)) .LE. reach ) THEN
      Eul_ind_1u(1) = i
      DO j = Eul_ind_1u(1), N11
        IF ( (x1u(j) - yb(Lag_ind,1)) .GT. reach ) THEN
          Eul_ind_1u(2) = (j-1)
          EXIT FIND_1u
        ELSE IF ( j == N11 ) THEN
          Eul_ind_1u(2) = j
          EXIT FIND_1u
        END IF
      END DO
    ELSE
      Eul_ind_1u = -1
    END IF
  END DO FIND_1u

  !--- compare the Lagrangian coordinate's second entry against the x2v grid (for vel 2) ---
  FIND_2v: DO i = S22, N22 !S22B, N22B?
    !--- add the first and last Eulerian coordinate's index that is within 'reach' to the list ---
    IF ( abs(yb(Lag_ind,2) - x2v(i)) .LE. reach ) THEN
      Eul_ind_2v(1) = i
      DO j = Eul_ind_2v(1), N22
        IF ( (x2v(j) - yb(Lag_ind,2)) .GT. reach ) THEN
          Eul_ind_2v(2) = (j-1)
          EXIT FIND_2v
        ELSE IF ( j == N22 ) THEN
          Eul_ind_2v(2) = j
          EXIT FIND_2v
        END IF
      END DO
    ELSE
      Eul_ind_2v = -1
    END IF
  END DO FIND_2v
  
  IF (dimens == 3) THEN
  !--- comppare the Lagrangian coordinate's third entry against the x3w grid (for vel 3) ---
  FIND_3w: DO i = S33, N33 !S33B, N33B?
    !--- add the first and last Eulerian coordinate's index that is within 'reach' to the list ---
    IF ( abs(yb(Lag_ind,3) - x3w(i)) .LE. reach ) THEN
      Eul_ind_3w(1) = i
      DO j = Eul_ind_3w(1), N33
        IF ( (x3w(j) - yb(Lag_ind,3)) .GT. reach ) THEN
          Eul_ind_3w(2) = (j-1)
          EXIT FIND_3w
        ELSE IF ( j == N33 ) THEN
          Eul_ind_3w(2) = j
          EXIT FIND_3w
        END IF
      END DO
    ELSE
      Eul_ind_3w = -1
    END IF
  END DO FIND_3w
  END IF

  
 


  END SUBROUTINE find_vicinity_coordinates





  SUBROUTINE compute_distances_to_vicinity_coordinates(Lag_ind, Eul_ind_1p, Eul_ind_1u, &
                                                                Eul_ind_2p, Eul_ind_2v, &
                                                                Eul_ind_3p, Eul_ind_3w, &
                                                                dist_1p, dist_1u, &
                                                                dist_2p, dist_2v, &
                                                                dist_3p, dist_3w)
  IMPLICIT NONE

  INTEGER, INTENT(IN)               :: Lag_ind
  INTEGER, INTENT(IN), DIMENSION(2) :: Eul_ind_1p, Eul_ind_1u
  INTEGER, INTENT(IN), DIMENSION(2) :: Eul_ind_2p, Eul_ind_2v
  INTEGER, INTENT(IN), DIMENSION(2) :: Eul_ind_3p, Eul_ind_3w

  REAL,    INTENT(OUT),DIMENSION(Eul_ind_1p(2) - Eul_ind_1p(1) + 1)           :: dist_1p
  REAL,    INTENT(OUT),DIMENSION(Eul_ind_2p(2) - Eul_ind_2p(1) + 1)           :: dist_2p
  REAL,    INTENT(OUT),DIMENSION(Eul_ind_3p(2) - Eul_ind_3p(1) + 1), OPTIONAL :: dist_3p
  REAL,    INTENT(OUT),DIMENSION(Eul_ind_1u(2) - Eul_ind_1u(1) + 1)           :: dist_1u
  REAL,    INTENT(OUT),DIMENSION(Eul_ind_2v(2) - Eul_ind_2v(1) + 1)           :: dist_2v
  REAL,    INTENT(OUT),DIMENSION(Eul_ind_3w(2) - Eul_ind_3w(1) + 1), OPTIONAL :: dist_3w

  INTEGER                           :: i
  
  
  !write(*,*) eul_ind_1p(1),',',eul_ind_1p(2),',',eul_ind_2p(1),',',eul_ind_2p(2),',',eul_ind_3p(1),',',eul_ind_3p(2),',',&
  !       eul_ind_1u(1),',',eul_ind_1u(2),',',eul_ind_2v(1),',',eul_ind_2v(2),',',eul_ind_3w(1),',',eul_ind_3w(2)

  !============================================================================================================
  IF (Eul_ind_1p(1) .NE. -1) THEN
    DO i = Eul_ind_1p(1), Eul_ind_1p(2)
      dist_1p((i - Eul_ind_1p(1)) + 1) = yb(Lag_ind,1) - x1p(i) 
    END DO
  ELSE
    dist_1p = L1*L2*L3
  END IF
  !============================================================================================================
  IF (Eul_ind_2p(1) .NE. -1) THEN
    DO i = Eul_ind_2p(1), Eul_ind_2p(2)
      dist_2p((i - Eul_ind_2p(1)) + 1) = yb(Lag_ind,2) - x2p(i) 
    END DO
  ELSE
    dist_2p = L1*L2*L3
  END IF
  !============================================================================================================
  IF (dimens == 3 .AND. Eul_ind_3p(1) .NE. -1) THEN
    DO i = Eul_ind_3p(1), Eul_ind_3p(2)
      dist_3p((i - Eul_ind_3p(1)) + 1) = yb(Lag_ind,3) - x3p(i) 
    END DO
  ELSE IF (dimens == 3) THEN
    dist_3p = L1*L2*L3
  END IF
  !============================================================================================================
  IF (Eul_ind_1u(1) .NE. -1) THEN
    DO i = Eul_ind_1u(1), Eul_ind_1u(2)
      dist_1u((i - Eul_ind_1u(1)) + 1) = yb(Lag_ind,1) - x1u(i) 
    END DO
  ELSE
    dist_1u = L1*L2*L3
  END IF
  !============================================================================================================
  IF ( Eul_ind_2v(1) .NE. -1) THEN
    DO i = Eul_ind_2v(1), Eul_ind_2v(2)
      dist_2v((i - Eul_ind_2v(1)) + 1) = yb(Lag_ind,2) - x2v(i) 
    END DO
  ELSE
    dist_2v = L1*L2*L3
  END IF
  !============================================================================================================
  IF (dimens == 3 .AND. Eul_ind_3w(1) .NE. -1) THEN
    DO i = Eul_ind_3w(1), Eul_ind_3w(2)
      dist_3w((i - Eul_ind_3w(1)) + 1) = yb(Lag_ind,3) - x3w(i) 
    END DO
  ELSE IF (dimens == 3) THEN
    dist_3w = L1*L2*L3
  END IF

  
 


  END SUBROUTINE compute_distances_to_vicinity_coordinates





  SUBROUTINE calculate_node_volumes(elem, elem_vol, f1, f2, f3, g1, g2, g3, h1, h2, h3, x1, x2, x3, y1, y2, y3)
  
    IMPLICIT NONE

    INTEGER                     :: node1, node2, node3
    REAL                        :: xx1, xx2, xx3, yy1, yy2, yy3
    REAL                        :: ff1, ff2, ff3, gg1, gg2, gg3, hh1, hh2, hh3
    REAL                        :: A
    
    INTEGER,INTENT(IN)          :: elem    
    REAL,INTENT(OUT),OPTIONAL   :: elem_vol
    REAL,INTENT(OUT),OPTIONAL   :: f1, f2, f3, g1, g2, g3, h1, h2, h3, x1, x2, x3, y1, y2, y3

      node1 = elems(elem,1)
      node2 = elems(elem,2)
      node3 = elems(elem,3)
    
      xx1 = yb( node1 , 1 )
      yy1 = yb( node1 , 2 )
      xx2 = yb( node2 , 1 )
      yy2 = yb( node2 , 2 )
      xx3 = yb( node3 , 1 )
      yy3 = yb( node3 , 2 )

      !write(*,*) i, elems(i,1),elems(i,2),elems(i,3)
      !write(*,*) x1,y1,x2,y2,x3,y3

      ff1 = xx2*yy3 - xx3*yy2
      ff2 = xx3*yy1 - xx1*yy3
      ff3 = xx1*yy2 - xx2*yy1

      gg1 = yy2 - yy3
      gg2 = yy3 - yy1
      gg3 = yy1 - yy2

      hh1 = xx3 - xx2
      hh2 = xx1 - xx3
      hh3 = xx2 - xx1

      A   = (ff1 + ff2 + ff3)/2.

    IF (PRESENT(elem_vol) ) THEN
      elem_vol  = A
      IF (PRESENT(f1) ) THEN
        f1 = ff1
        f2 = ff2
        f3 = ff3
        g1 = gg1
        g2 = gg2 
        g3 = gg3
        h1 = hh1
        h2 = hh2
        h3 = hh3
     END IF
   END IF
   IF (PRESENT(x1) ) THEN
      x1 = xx1
      x2 = xx2
      x3 = xx3
      y1 = yy1
      y2 = yy2
      y3 = yy3
    END IF

    node_vol(node1) = node_vol(node1) + A/3.
    node_vol(node2) = node_vol(node2) + A/3.
    node_vol(node3) = node_vol(node3) + A/3.


  END SUBROUTINE calculate_node_volumes






  SUBROUTINE update_boundary
    
    IMPLICIT NONE
 
    INTEGER          :: i,j           
    !CHARACTER(LEN=5) :: timestep_char !*** debugging

    !*** debugging
    !IF (substep == 1 .AND. (MOD(timestep,10) == 0)) then
    !CALL num_to_string(5,(timestep/10),timestep_char)
    !open(unit=54,file=('bound_coords_'//timestep_char//'.txt'),status='unknown')
    !do i=1,M_bound
    !write(54,*), (yb(i,j), j=1,2)
    !end do    
    !endif
   
    !ub = ub*(mu_fluid*Re/(rho_fluid*U_ref)) ! bbecsek 150507: make dimensional

    !$omp parallel do collapse(2)
    DO j = 1, dimens
      UPD_LOOP: DO i = 1, M_bound
        IF ( ANY( i == NINT(ebcs(:,1) ) ) ) cycle UPD_LOOP
        yb(i,j) = yb(i,j) + dtime*aRK(substep)*ub(i,j) + dtime*bRK(substep)*container(i,j) 
      END DO UPD_LOOP
    END DO
    !$omp end parallel do


    !--- broadcast boundary coordinate across all processes ---
    !CALL MPI_IBCAST(yb,SIZE(yb),MPI_REAL8,0,COMM_CART,req_bcast,merror)
    !CALL MPI_WAIT(req_bcast,status,merror)
  
  END SUBROUTINE update_boundary



  !< subrouine that computes displacement field, called in mod_timeint in each sub-timestep
  SUBROUTINE calculate_displacements
    
    IMPLICIT NONE

    INTEGER          :: i,j
    
    db = 0.

    !$omp parallel do collapse(2)
    DO j = 1, dimens
      DISP_LOOP: DO i = 1, M_bound
        IF ( ANY( i == NINT(ebcs(:,1) ) ) ) cycle DISP_LOOP
        db(i,j) = yb(i,j) - xb(i,j)
      END DO DISP_LOOP
    END DO
    !$omp end parallel do

    !--- broadcast boundary displacement across all processes ---
    !CALL MPI_IBCAST(db,SIZE(db),MPI_REAL8,0,COMM_CART,req_bcast,merror)
    !CALL MPI_WAIT(req_bcast,status,merror)
  
  END SUBROUTINE calculate_displacements


  !< subroutine to set up matrix B to link nodal DOFs and strain, calculate strain/stress and set up global stiffness matrix
  SUBROUTINE fe_setup_triangular_3node
  
    IMPLICIT NONE

    CALL setup_stiffness_matrix
    CALL setup_mass_matrix
    !CALL setup_damping_matrix

  END SUBROUTINE fe_setup_triangular_3node


  SUBROUTINE compute_strain_stress

    IMPLICIT NONE

    REAL                 :: f1,f2,f3
    REAL                 :: g1,g2,g3
    REAL                 :: h1,h2,h3
    REAL                 :: A
    REAL, DIMENSION(3,6) :: B_transp

    INTEGER              :: i

    strain        = 0. 
    stress        = 0.
    node_vol      = 0.



    DO i = 1, M_elems

      CALL calculate_node_volumes(i, A, f1, f2, f3, g1, g2, g3, h1, h2, h3)
      !write(*,*) node_vol


      B_transp = 0.
      !--- set up matrix B^T that links nodal DOFs with the strain ---
      B_transp(1,:) = (/g1,0.,g2,0.,g3,0./)
      B_transp(2,:) = (/0.,h1,0.,h2,0.,h3/)
      B_transp(3,:) = (/h1,g1,h2,g2,h3,g3/)

      B_transp = (1./(2.*A))*B_transp


      !--- calculate strain of i-th element (B^T*d) ---
      strain(i,1) = B_transp(1,1)*db(elems(i,1),1) + B_transp(1,2)*db(elems(i,1),2) &
                  + B_transp(1,3)*db(elems(i,2),1) + B_transp(1,4)*db(elems(i,2),2) &
                  + B_transp(1,5)*db(elems(i,3),1) + B_transp(1,6)*db(elems(i,3),2)
      strain(i,2) = B_transp(2,1)*db(elems(i,1),1) + B_transp(2,2)*db(elems(i,1),2) &
                  + B_transp(2,3)*db(elems(i,2),1) + B_transp(2,4)*db(elems(i,2),2) &
                  + B_transp(2,5)*db(elems(i,3),1) + B_transp(2,6)*db(elems(i,3),2)
      strain(i,3) = B_transp(3,1)*db(elems(i,1),1) + B_transp(3,2)*db(elems(i,1),2) &
                  + B_transp(3,3)*db(elems(i,2),1) + B_transp(3,4)*db(elems(i,2),2) &
                  + B_transp(3,5)*db(elems(i,3),1) + B_transp(3,6)*db(elems(i,3),2)

      !--- calculate stress of i-th element (C*strain) ---
      stress(i,1) = C(1,1)*strain(i,1) + C(1,2)*strain(i,2) + C(1,3)*strain(i,3) 
      stress(i,2) = C(2,1)*strain(i,1) + C(2,2)*strain(i,2) + C(2,3)*strain(i,3) 
      stress(i,3) = C(3,1)*strain(i,1) + C(3,2)*strain(i,2) + C(3,3)*strain(i,3) 


    END DO
 
  END SUBROUTINE compute_strain_stress



  SUBROUTINE setup_stiffness_matrix

    IMPLICIT NONE

    !REAL                 :: x1,x2,x3,y1,y2,y3
    REAL                 :: f1,f2,f3
    REAL                 :: g1,g2,g3
    REAL                 :: h1,h2,h3
    REAL                 :: A
    REAL, DIMENSION(3,6) :: B_transp
    REAL, DIMENSION(6,6) :: k_stiff_local
    REAL                 :: k_forsend, k_reduced

    INTEGER          :: glob_dof_x1, glob_dof_x2, glob_dof_x3, glob_dof_y1, glob_dof_y2, glob_dof_y3 
    INTEGER          :: i,j


    k_stiff       = 0.
    strain        = 0. 
    stress        = 0.
    node_vol      = 0.

    !$omp parallel do private(B_transp,A,f1,f2,f3,g1,g2,g3,h1,h2,h3,k_stiff_local,glob_dof_x1,glob_dof_x2,glob_dof_x3,glob_dof_y1,glob_dof_y2,glob_dof_y3) reduction(+:k_stiff)
    DO i = 1,M_elems


      CALL calculate_node_volumes(i, A, f1, f2, f3, g1, g2, g3, h1, h2, h3)
      !write(*,*) node_vol


      B_transp = 0.
      !--- set up matrix B^T that links nodal DOFs with the strain ---
      B_transp(1,:) = (/g1,0.,g2,0.,g3,0./)
      B_transp(2,:) = (/0.,h1,0.,h2,0.,h3/)
      B_transp(3,:) = (/h1,g1,h2,g2,h3,g3/)

      B_transp = (1./(2.*A))*B_transp


      !--- calculate strain of i-th element (B^T*d) ---
      strain(i,1) = B_transp(1,1)*db(elems(i,1),1) + B_transp(1,2)*db(elems(i,1),2) &
                  + B_transp(1,3)*db(elems(i,2),1) + B_transp(1,4)*db(elems(i,2),2) &
                  + B_transp(1,5)*db(elems(i,3),1) + B_transp(1,6)*db(elems(i,3),2)
      strain(i,2) = B_transp(2,1)*db(elems(i,1),1) + B_transp(2,2)*db(elems(i,1),2) &
                  + B_transp(2,3)*db(elems(i,2),1) + B_transp(2,4)*db(elems(i,2),2) &
                  + B_transp(2,5)*db(elems(i,3),1) + B_transp(2,6)*db(elems(i,3),2)
      strain(i,3) = B_transp(3,1)*db(elems(i,1),1) + B_transp(3,2)*db(elems(i,1),2) &
                  + B_transp(3,3)*db(elems(i,2),1) + B_transp(3,4)*db(elems(i,2),2) &
                  + B_transp(3,5)*db(elems(i,3),1) + B_transp(3,6)*db(elems(i,3),2)

      !--- calculate stress of i-th element (C*strain) ---
      stress(i,1) = C(1,1)*strain(i,1) + C(1,2)*strain(i,2) + C(1,3)*strain(i,3) 
      stress(i,2) = C(2,1)*strain(i,1) + C(2,2)*strain(i,2) + C(2,3)*strain(i,3) 
      stress(i,3) = C(3,1)*strain(i,1) + C(3,2)*strain(i,2) + C(3,3)*strain(i,3) 
      
      k_stiff_local = 0.

      k_stiff_local = MATMUL(TRANSPOSE(B_transp),MATMUL(C,B_transp))
      !--- calculate local stiffness matrix ---
      !--- x-component of node 1 ---
!      k_stiff_local(1,1) =   (B_transp(1,1)*( C(1,1)*B_transp(1,1) + C(1,2)*B_transp(2,1) + C(1,3)*B_transp(3,1) ) &
!                            + B_transp(2,1)*( C(2,1)*B_transp(1,1) + C(2,2)*B_transp(2,1) + C(2,3)*B_transp(3,1) ) &
!                            + B_transp(3,1)*( C(3,1)*B_transp(1,1) + C(3,2)*B_transp(2,1) + C(3,3)*B_transp(3,1) ) )
!      
!      k_stiff_local(1,2) =   (B_transp(1,1)*( C(1,1)*B_transp(1,2) + C(1,2)*B_transp(2,2) + C(1,3)*B_transp(3,2) ) &
!                            + B_transp(2,1)*( C(2,1)*B_transp(1,2) + C(2,2)*B_transp(2,2) + C(2,3)*B_transp(3,2) ) &
!                            + B_transp(3,1)*( C(3,1)*B_transp(1,2) + C(3,2)*B_transp(2,2) + C(3,3)*B_transp(3,2) ) )
! 
!      k_stiff_local(1,3) =   (B_transp(1,1)*( C(1,1)*B_transp(1,3) + C(1,2)*B_transp(2,3) + C(1,3)*B_transp(3,3) ) &
!                            + B_transp(2,1)*( C(2,1)*B_transp(1,3) + C(2,2)*B_transp(2,3) + C(2,3)*B_transp(3,3) ) &
!                            + B_transp(3,1)*( C(3,1)*B_transp(1,3) + C(3,2)*B_transp(2,3) + C(3,3)*B_transp(3,3) ) )
! 
!      k_stiff_local(1,4) =   (B_transp(1,1)*( C(1,1)*B_transp(1,4) + C(1,2)*B_transp(2,4) + C(1,3)*B_transp(3,4) ) &
!                            + B_transp(2,1)*( C(2,1)*B_transp(1,4) + C(2,2)*B_transp(2,4) + C(2,3)*B_transp(3,4) ) &
!                            + B_transp(3,1)*( C(3,1)*B_transp(1,4) + C(3,2)*B_transp(2,4) + C(3,3)*B_transp(3,4) ) )
! 
!      k_stiff_local(1,5) =   (B_transp(1,1)*( C(1,1)*B_transp(1,5) + C(1,2)*B_transp(2,5) + C(1,3)*B_transp(3,5) ) &
!                            + B_transp(2,1)*( C(2,1)*B_transp(1,5) + C(2,2)*B_transp(2,5) + C(2,3)*B_transp(3,5) ) &
!                            + B_transp(3,1)*( C(3,1)*B_transp(1,5) + C(3,2)*B_transp(2,5) + C(3,3)*B_transp(3,5) ) )
! 
!      k_stiff_local(1,6) =   (B_transp(1,1)*( C(1,1)*B_transp(1,6) + C(1,2)*B_transp(2,6) + C(1,3)*B_transp(3,6) ) &
!                            + B_transp(2,1)*( C(2,1)*B_transp(1,6) + C(2,2)*B_transp(2,6) + C(2,3)*B_transp(3,6) ) &
!                            + B_transp(3,1)*( C(3,1)*B_transp(1,6) + C(3,2)*B_transp(2,6) + C(3,3)*B_transp(3,6) ) )
! 
!      !--- y-component of node 1 ---
!      k_stiff_local(2,1) =   (B_transp(1,2)*( C(1,1)*B_transp(1,1) + C(1,2)*B_transp(2,1) + C(1,3)*B_transp(3,1) ) &
!                            + B_transp(2,2)*( C(2,1)*B_transp(1,1) + C(2,2)*B_transp(2,1) + C(2,3)*B_transp(3,1) ) &
!                            + B_transp(3,2)*( C(3,1)*B_transp(1,1) + C(3,2)*B_transp(2,1) + C(3,3)*B_transp(3,1) ) )
! 
!      k_stiff_local(2,2) =   (B_transp(1,2)*( C(1,1)*B_transp(1,2) + C(1,2)*B_transp(2,2) + C(1,3)*B_transp(3,2) ) &
!                            + B_transp(2,2)*( C(2,1)*B_transp(1,2) + C(2,2)*B_transp(2,2) + C(2,3)*B_transp(3,2) ) &
!                            + B_transp(3,2)*( C(3,1)*B_transp(1,2) + C(3,2)*B_transp(2,2) + C(3,3)*B_transp(3,2) ) )
! 
!      k_stiff_local(2,3) =   (B_transp(1,2)*( C(1,1)*B_transp(1,3) + C(1,2)*B_transp(2,3) + C(1,3)*B_transp(3,3) ) &
!                            + B_transp(2,2)*( C(2,1)*B_transp(1,3) + C(2,2)*B_transp(2,3) + C(2,3)*B_transp(3,3) ) &
!                            + B_transp(3,2)*( C(3,1)*B_transp(1,3) + C(3,2)*B_transp(2,3) + C(3,3)*B_transp(3,3) ) )
! 
!      k_stiff_local(2,4) =   (B_transp(1,2)*( C(1,1)*B_transp(1,4) + C(1,2)*B_transp(2,4) + C(1,3)*B_transp(3,4) ) &
!                            + B_transp(2,2)*( C(2,1)*B_transp(1,4) + C(2,2)*B_transp(2,4) + C(2,3)*B_transp(3,4) ) &
!                            + B_transp(3,2)*( C(3,1)*B_transp(1,4) + C(3,2)*B_transp(2,4) + C(3,3)*B_transp(3,4) ) )
! 
!      k_stiff_local(2,5) =   (B_transp(1,2)*( C(1,1)*B_transp(1,5) + C(1,2)*B_transp(2,5) + C(1,3)*B_transp(3,5) ) &
!                            + B_transp(2,2)*( C(2,1)*B_transp(1,5) + C(2,2)*B_transp(2,5) + C(2,3)*B_transp(3,5) ) &
!                            + B_transp(3,2)*( C(3,1)*B_transp(1,5) + C(3,2)*B_transp(2,5) + C(3,3)*B_transp(3,5) ) )
! 
!      k_stiff_local(2,6) =   (B_transp(1,2)*( C(1,1)*B_transp(1,6) + C(1,2)*B_transp(2,6) + C(1,3)*B_transp(3,6) ) &
!                            + B_transp(2,2)*( C(2,1)*B_transp(1,6) + C(2,2)*B_transp(2,6) + C(2,3)*B_transp(3,6) ) &
!                            + B_transp(3,2)*( C(3,1)*B_transp(1,6) + C(3,2)*B_transp(2,6) + C(3,3)*B_transp(3,6) ) )
! 
!      !--- x-component of node 2 ---
!      k_stiff_local(3,1) =   (B_transp(1,3)*( C(1,1)*B_transp(1,1) + C(1,2)*B_transp(2,1) + C(1,3)*B_transp(3,1) ) &
!                            + B_transp(2,3)*( C(2,1)*B_transp(1,1) + C(2,2)*B_transp(2,1) + C(2,3)*B_transp(3,1) ) &
!                            + B_transp(3,3)*( C(3,1)*B_transp(1,1) + C(3,2)*B_transp(2,1) + C(3,3)*B_transp(3,1) ) )
! 
!      k_stiff_local(3,2) =   (B_transp(1,3)*( C(1,1)*B_transp(1,2) + C(1,2)*B_transp(2,2) + C(1,3)*B_transp(3,2) ) &
!                            + B_transp(2,3)*( C(2,1)*B_transp(1,2) + C(2,2)*B_transp(2,2) + C(2,3)*B_transp(3,2) ) &
!                            + B_transp(3,3)*( C(3,1)*B_transp(1,2) + C(3,2)*B_transp(2,2) + C(3,3)*B_transp(3,2) ) )
! 
!      k_stiff_local(3,3) =   (B_transp(1,3)*( C(1,1)*B_transp(1,3) + C(1,2)*B_transp(2,3) + C(1,3)*B_transp(3,3) ) &
!                            + B_transp(2,3)*( C(2,1)*B_transp(1,3) + C(2,2)*B_transp(2,3) + C(2,3)*B_transp(3,3) ) &
!                            + B_transp(3,3)*( C(3,1)*B_transp(1,3) + C(3,2)*B_transp(2,3) + C(3,3)*B_transp(3,3) ) )
! 
!      k_stiff_local(3,4) =   (B_transp(1,3)*( C(1,1)*B_transp(1,4) + C(1,2)*B_transp(2,4) + C(1,3)*B_transp(3,4) ) &
!                            + B_transp(2,3)*( C(2,1)*B_transp(1,4) + C(2,2)*B_transp(2,4) + C(2,3)*B_transp(3,4) ) &
!                            + B_transp(3,3)*( C(3,1)*B_transp(1,4) + C(3,2)*B_transp(2,4) + C(3,3)*B_transp(3,4) ) )
! 
!      k_stiff_local(3,5) =   (B_transp(1,3)*( C(1,1)*B_transp(1,5) + C(1,2)*B_transp(2,5) + C(1,3)*B_transp(3,5) ) &
!                            + B_transp(2,3)*( C(2,1)*B_transp(1,5) + C(2,2)*B_transp(2,5) + C(2,3)*B_transp(3,5) ) &
!                            + B_transp(3,3)*( C(3,1)*B_transp(1,5) + C(3,2)*B_transp(2,5) + C(3,3)*B_transp(3,5) ) )
! 
!      k_stiff_local(3,6) =   (B_transp(1,3)*( C(1,1)*B_transp(1,6) + C(1,2)*B_transp(2,6) + C(1,3)*B_transp(3,6) ) &
!                            + B_transp(2,3)*( C(2,1)*B_transp(1,6) + C(2,2)*B_transp(2,6) + C(2,3)*B_transp(3,6) ) &
!                            + B_transp(3,3)*( C(3,1)*B_transp(1,6) + C(3,2)*B_transp(2,6) + C(3,3)*B_transp(3,6) ) )
! 
!      !--- y-component of node 2 ---
!      k_stiff_local(4,1) =   (B_transp(1,4)*( C(1,1)*B_transp(1,1) + C(1,2)*B_transp(2,1) + C(1,3)*B_transp(3,1) ) &
!                            + B_transp(2,4)*( C(2,1)*B_transp(1,1) + C(2,2)*B_transp(2,1) + C(2,3)*B_transp(3,1) ) &
!                            + B_transp(3,4)*( C(3,1)*B_transp(1,1) + C(3,2)*B_transp(2,1) + C(3,3)*B_transp(3,1) ) )
! 
!      k_stiff_local(4,2) =   (B_transp(1,4)*( C(1,1)*B_transp(1,2) + C(1,2)*B_transp(2,2) + C(1,3)*B_transp(3,2) ) &
!                            + B_transp(2,4)*( C(2,1)*B_transp(1,2) + C(2,2)*B_transp(2,2) + C(2,3)*B_transp(3,2) ) &
!                            + B_transp(3,4)*( C(3,1)*B_transp(1,2) + C(3,2)*B_transp(2,2) + C(3,3)*B_transp(3,2) ) )
! 
!      k_stiff_local(4,3) =   (B_transp(1,4)*( C(1,1)*B_transp(1,3) + C(1,2)*B_transp(2,3) + C(1,3)*B_transp(3,3) ) &
!                            + B_transp(2,4)*( C(2,1)*B_transp(1,3) + C(2,2)*B_transp(2,3) + C(2,3)*B_transp(3,3) ) &
!                            + B_transp(3,4)*( C(3,1)*B_transp(1,3) + C(3,2)*B_transp(2,3) + C(3,3)*B_transp(3,3) ) )
! 
!      k_stiff_local(4,4) =   (B_transp(1,4)*( C(1,1)*B_transp(1,4) + C(1,2)*B_transp(2,4) + C(1,3)*B_transp(3,4) ) &
!                            + B_transp(2,4)*( C(2,1)*B_transp(1,4) + C(2,2)*B_transp(2,4) + C(2,3)*B_transp(3,4) ) &
!                            + B_transp(3,4)*( C(3,1)*B_transp(1,4) + C(3,2)*B_transp(2,4) + C(3,3)*B_transp(3,4) ) )
! 
!      k_stiff_local(4,5) =   (B_transp(1,4)*( C(1,1)*B_transp(1,5) + C(1,2)*B_transp(2,5) + C(1,3)*B_transp(3,5) ) &
!                            + B_transp(2,4)*( C(2,1)*B_transp(1,5) + C(2,2)*B_transp(2,5) + C(2,3)*B_transp(3,5) ) &
!                            + B_transp(3,4)*( C(3,1)*B_transp(1,5) + C(3,2)*B_transp(2,5) + C(3,3)*B_transp(3,5) ) )
! 
!      k_stiff_local(4,6) =   (B_transp(1,4)*( C(1,1)*B_transp(1,6) + C(1,2)*B_transp(2,6) + C(1,3)*B_transp(3,6) ) &
!                            + B_transp(2,4)*( C(2,1)*B_transp(1,6) + C(2,2)*B_transp(2,6) + C(2,3)*B_transp(3,6) ) &
!                            + B_transp(3,4)*( C(3,1)*B_transp(1,6) + C(3,2)*B_transp(2,6) + C(3,3)*B_transp(3,6) ) )
! 
!      !--- x-component of node 3 ---
!      k_stiff_local(5,1) =   (B_transp(1,5)*( C(1,1)*B_transp(1,1) + C(1,2)*B_transp(2,1) + C(1,3)*B_transp(3,1) ) &
!                            + B_transp(2,5)*( C(2,1)*B_transp(1,1) + C(2,2)*B_transp(2,1) + C(2,3)*B_transp(3,1) ) &
!                            + B_transp(3,5)*( C(3,1)*B_transp(1,1) + C(3,2)*B_transp(2,1) + C(3,3)*B_transp(3,1) ) )
! 
!      k_stiff_local(5,2) =   (B_transp(1,5)*( C(1,1)*B_transp(1,2) + C(1,2)*B_transp(2,2) + C(1,3)*B_transp(3,2) ) &
!                            + B_transp(2,5)*( C(2,1)*B_transp(1,2) + C(2,2)*B_transp(2,2) + C(2,3)*B_transp(3,2) ) &
!                            + B_transp(3,5)*( C(3,1)*B_transp(1,2) + C(3,2)*B_transp(2,2) + C(3,3)*B_transp(3,2) ) )
! 
!      k_stiff_local(5,3) =   (B_transp(1,5)*( C(1,1)*B_transp(1,3) + C(1,2)*B_transp(2,3) + C(1,3)*B_transp(3,3) ) &
!                            + B_transp(2,5)*( C(2,1)*B_transp(1,3) + C(2,2)*B_transp(2,3) + C(2,3)*B_transp(3,3) ) &
!                            + B_transp(3,5)*( C(3,1)*B_transp(1,3) + C(3,2)*B_transp(2,3) + C(3,3)*B_transp(3,3) ) )
! 
!      k_stiff_local(5,4) =   (B_transp(1,5)*( C(1,1)*B_transp(1,4) + C(1,2)*B_transp(2,4) + C(1,3)*B_transp(3,4) ) &
!                            + B_transp(2,5)*( C(2,1)*B_transp(1,4) + C(2,2)*B_transp(2,4) + C(2,3)*B_transp(3,4) ) &
!                            + B_transp(3,5)*( C(3,1)*B_transp(1,4) + C(3,2)*B_transp(2,4) + C(3,3)*B_transp(3,4) ) ) 
!     
!      k_stiff_local(5,5) =   (B_transp(1,5)*( C(1,1)*B_transp(1,5) + C(1,2)*B_transp(2,5) + C(1,3)*B_transp(3,5) ) &
!                            + B_transp(2,5)*( C(2,1)*B_transp(1,5) + C(2,2)*B_transp(2,5) + C(2,3)*B_transp(3,5) ) &
!                            + B_transp(3,5)*( C(3,1)*B_transp(1,5) + C(3,2)*B_transp(2,5) + C(3,3)*B_transp(3,5) ) )
! 
!      k_stiff_local(5,6) =   (B_transp(1,5)*( C(1,1)*B_transp(1,6) + C(1,2)*B_transp(2,6) + C(1,3)*B_transp(3,6) ) &
!                            + B_transp(2,5)*( C(2,1)*B_transp(1,6) + C(2,2)*B_transp(2,6) + C(2,3)*B_transp(3,6) ) &
!                            + B_transp(3,5)*( C(3,1)*B_transp(1,6) + C(3,2)*B_transp(2,6) + C(3,3)*B_transp(3,6) ) )
! 
!      !--- y-component of node 3 ---
!      k_stiff_local(6,1) =   (B_transp(1,6)*( C(1,1)*B_transp(1,1) + C(1,2)*B_transp(2,1) + C(1,3)*B_transp(3,1) ) &
!                            + B_transp(2,6)*( C(2,1)*B_transp(1,1) + C(2,2)*B_transp(2,1) + C(2,3)*B_transp(3,1) ) &
!                            + B_transp(3,6)*( C(3,1)*B_transp(1,1) + C(3,2)*B_transp(2,1) + C(3,3)*B_transp(3,1) ) )
! 
!      k_stiff_local(6,2) =   (B_transp(1,6)*( C(1,1)*B_transp(1,2) + C(1,2)*B_transp(2,2) + C(1,3)*B_transp(3,2) ) &
!                            + B_transp(2,6)*( C(2,1)*B_transp(1,2) + C(2,2)*B_transp(2,2) + C(2,3)*B_transp(3,2) ) &
!                            + B_transp(3,6)*( C(3,1)*B_transp(1,2) + C(3,2)*B_transp(2,2) + C(3,3)*B_transp(3,2) ) )
! 
!      k_stiff_local(6,3) =   (B_transp(1,6)*( C(1,1)*B_transp(1,3) + C(1,2)*B_transp(2,3) + C(1,3)*B_transp(3,3) ) &
!                            + B_transp(2,6)*( C(2,1)*B_transp(1,3) + C(2,2)*B_transp(2,3) + C(2,3)*B_transp(3,3) ) &
!                            + B_transp(3,6)*( C(3,1)*B_transp(1,3) + C(3,2)*B_transp(2,3) + C(3,3)*B_transp(3,3) ) )
! 
!      k_stiff_local(6,4) =   (B_transp(1,6)*( C(1,1)*B_transp(1,4) + C(1,2)*B_transp(2,4) + C(1,3)*B_transp(3,4) ) &
!                            + B_transp(2,6)*( C(2,1)*B_transp(1,4) + C(2,2)*B_transp(2,4) + C(2,3)*B_transp(3,4) ) &
!                            + B_transp(3,6)*( C(3,1)*B_transp(1,4) + C(3,2)*B_transp(2,4) + C(3,3)*B_transp(3,4) ) )
! 
!      k_stiff_local(6,5) =   (B_transp(1,6)*( C(1,1)*B_transp(1,5) + C(1,2)*B_transp(2,5) + C(1,3)*B_transp(3,5) ) &
!                            + B_transp(2,6)*( C(2,1)*B_transp(1,5) + C(2,2)*B_transp(2,5) + C(2,3)*B_transp(3,5) ) &
!                            + B_transp(3,6)*( C(3,1)*B_transp(1,5) + C(3,2)*B_transp(2,5) + C(3,3)*B_transp(3,5) ) )
! 
!      k_stiff_local(6,6) =   (B_transp(1,6)*( C(1,1)*B_transp(1,6) + C(1,2)*B_transp(2,6) + C(1,3)*B_transp(3,6) ) &
!                            + B_transp(2,6)*( C(2,1)*B_transp(1,6) + C(2,2)*B_transp(2,6) + C(2,3)*B_transp(3,6) ) &
!                            + B_transp(3,6)*( C(3,1)*B_transp(1,6) + C(3,2)*B_transp(2,6) + C(3,3)*B_transp(3,6) ) )


      k_stiff_local = A*k_stiff_local

      !--- calculate global stiffness matrix ---
      !--- global DOFs ---
      glob_dof_x1 = elems(i,1)*2 - 1
      glob_dof_y1 = elems(i,1)*2
      glob_dof_x2 = elems(i,2)*2 - 1
      glob_dof_y2 = elems(i,2)*2
      glob_dof_x3 = elems(i,3)*2 - 1 
      glob_dof_y3 = elems(i,3)*2

        

      !--- x-component of node 1 ---
      k_stiff(glob_dof_x1, glob_dof_x1) = k_stiff(glob_dof_x1, glob_dof_x1) + k_stiff_local(1,1)
      
      k_stiff(glob_dof_x1, glob_dof_y1) = k_stiff(glob_dof_x1, glob_dof_y1) + k_stiff_local(1,2)
                                           
      k_stiff(glob_dof_x1, glob_dof_x2) = k_stiff(glob_dof_x1, glob_dof_x2) + k_stiff_local(1,3)
                                           
      k_stiff(glob_dof_x1, glob_dof_y2) = k_stiff(glob_dof_x1, glob_dof_y2) + k_stiff_local(1,4)
                                     
      k_stiff(glob_dof_x1, glob_dof_x3) = k_stiff(glob_dof_x1, glob_dof_x3) + k_stiff_local(1,5)
                                     
      k_stiff(glob_dof_x1, glob_dof_y3) = k_stiff(glob_dof_x1, glob_dof_y3) + k_stiff_local(1,6)
                                    
                                       
      !--- y-component of node 1 ---
      k_stiff(glob_dof_y1, glob_dof_x1) = k_stiff(glob_dof_y1, glob_dof_x1) + k_stiff_local(2,1)
                                          
      k_stiff(glob_dof_y1, glob_dof_y1) = k_stiff(glob_dof_y1, glob_dof_y1) + k_stiff_local(2,2)
                                           
      k_stiff(glob_dof_y1, glob_dof_x2) = k_stiff(glob_dof_y1, glob_dof_x2) + k_stiff_local(2,3)
                                           
      k_stiff(glob_dof_y1, glob_dof_y2) = k_stiff(glob_dof_y1, glob_dof_y2) + k_stiff_local(2,4)
                                     
      k_stiff(glob_dof_y1, glob_dof_x3) = k_stiff(glob_dof_y1, glob_dof_x3) + k_stiff_local(2,5)
                                     
      k_stiff(glob_dof_y1, glob_dof_y3) = k_stiff(glob_dof_y1, glob_dof_y3) + k_stiff_local(2,6)
                                     

      !--- x-component of node 2 ---
      k_stiff(glob_dof_x2, glob_dof_x1) = k_stiff(glob_dof_x2, glob_dof_x1) + k_stiff_local(3,1)
                                           
      k_stiff(glob_dof_x2, glob_dof_y1) = k_stiff(glob_dof_x2, glob_dof_y1) + k_stiff_local(3,2)
                                           
      k_stiff(glob_dof_x2, glob_dof_x2) = k_stiff(glob_dof_x2, glob_dof_x2) + k_stiff_local(3,3)
                                            
      k_stiff(glob_dof_x2, glob_dof_y2) = k_stiff(glob_dof_x2, glob_dof_y2) + k_stiff_local(3,4)
                                    
      k_stiff(glob_dof_x2, glob_dof_x3) = k_stiff(glob_dof_x2, glob_dof_x3) + k_stiff_local(3,5)
                                     
      k_stiff(glob_dof_x2, glob_dof_y3) = k_stiff(glob_dof_x2, glob_dof_y3) + k_stiff_local(3,6)
                                     

      !--- y-component of node 2 ---
      k_stiff(glob_dof_y2, glob_dof_x1) = k_stiff(glob_dof_y2, glob_dof_x1) + k_stiff_local(4,1)
                                           
      k_stiff(glob_dof_y2, glob_dof_y1) = k_stiff(glob_dof_y2, glob_dof_y1) + k_stiff_local(4,2) 
                                           
      k_stiff(glob_dof_y2, glob_dof_x2) = k_stiff(glob_dof_y2, glob_dof_x2) + k_stiff_local(4,3)
                                          
      k_stiff(glob_dof_y2, glob_dof_y2) = k_stiff(glob_dof_y2, glob_dof_y2) + k_stiff_local(4,4)
                                     
      k_stiff(glob_dof_y2, glob_dof_x3) = k_stiff(glob_dof_y2, glob_dof_x3) + k_stiff_local(4,5)
                                     
      k_stiff(glob_dof_y2, glob_dof_y3) = k_stiff(glob_dof_y2, glob_dof_y3) + k_stiff_local(4,6)
                                     

      !--- x-component of node 3 ---
      k_stiff(glob_dof_x3, glob_dof_x1) = k_stiff(glob_dof_x3, glob_dof_x1) + k_stiff_local(5,1)
                                           
      k_stiff(glob_dof_x3, glob_dof_y1) = k_stiff(glob_dof_x3, glob_dof_y1) + k_stiff_local(5,2)
                                           
      k_stiff(glob_dof_x3, glob_dof_x2) = k_stiff(glob_dof_x3, glob_dof_x2) + k_stiff_local(5,3)
                                           
      k_stiff(glob_dof_x3, glob_dof_y2) = k_stiff(glob_dof_x3, glob_dof_y2) + k_stiff_local(5,4)                

      k_stiff(glob_dof_x3, glob_dof_x3) = k_stiff(glob_dof_x3, glob_dof_x3) + k_stiff_local(5,5)
                                     
      k_stiff(glob_dof_x3, glob_dof_y3) = k_stiff(glob_dof_x3, glob_dof_y3) + k_stiff_local(5,6)
                                     

      !--- y-component of node 3 ---
      k_stiff(glob_dof_y3, glob_dof_x1) = k_stiff(glob_dof_y3, glob_dof_x1) + k_stiff_local(6,1)
                                           
      k_stiff(glob_dof_y3, glob_dof_y1) = k_stiff(glob_dof_y3, glob_dof_y1) + k_stiff_local(6,2)
                                          
      k_stiff(glob_dof_y3, glob_dof_x2) = k_stiff(glob_dof_y3, glob_dof_x2) + k_stiff_local(6,3)
                                           
      k_stiff(glob_dof_y3, glob_dof_y2) = k_stiff(glob_dof_y3, glob_dof_y2) + k_stiff_local(6,4)
                                     
      k_stiff(glob_dof_y3, glob_dof_x3) = k_stiff(glob_dof_y3, glob_dof_x3) + k_stiff_local(6,5)
                                     
      k_stiff(glob_dof_y3, glob_dof_y3) = k_stiff(glob_dof_y3, glob_dof_y3) + k_stiff_local(6,6)
                                     
    
    !write(*,*) (k_stiff_local(1,2).EQ.k_stiff_local(2,1))
    !write(*,*) (k_stiff_local(1,3).EQ.k_stiff_local(3,1)),(k_stiff_local(2,3).EQ.k_stiff_local(3,2))
    !write(*,*) (k_stiff_local(1,4).EQ.k_stiff_local(4,1)),(k_stiff_local(2,4).EQ.k_stiff_local(4,2)),(k_stiff_local(3,4).EQ.k_stiff_local(4,3))
    !write(*,*) k_stiff_local(1,4),k_stiff_local(4,1)
    !write(*,*) (k_stiff_local(1,5).EQ.k_stiff_local(5,1)),(k_stiff_local(2,5).EQ.k_stiff_local(5,2)),(k_stiff_local(3,5).EQ.k_stiff_local(5,3)),(k_stiff_local(4,5).EQ.k_stiff_local(5,4))
    !write(*,*) (k_stiff_local(1,6).EQ.k_stiff_local(6,1)),(k_stiff_local(2,6).EQ.k_stiff_local(6,2)),(k_stiff_local(3,6).EQ.k_stiff_local(6,3)),(k_stiff_local(4,6).EQ.k_stiff_local(6,4)),(k_stiff_local(5,6).EQ.k_stiff_local(6,5))
   

    END DO
    !$omp end parallel do

  END SUBROUTINE setup_stiffness_matrix




  SUBROUTINE setup_damping_matrix

    IMPLICIT NONE

    REAL     :: alpha = 0.2
    REAL     :: beta  = 0.8

    c_damp = alpha*m_mass + beta*k_stiff
    
  END SUBROUTINE setup_damping_matrix

 !> This is basically the same mass matrix as in setup_mass_matrix2 but requires less computation since the integrals for
 !! the triangular 3-node elements are known for general cases.
  SUBROUTINE setup_mass_matrix

    IMPLICIT NONE

    REAL                 :: A
    REAL, DIMENSION(6,6) :: m_mass_local

    INTEGER              :: glob_dof_x1, glob_dof_x2, glob_dof_x3, glob_dof_y1, glob_dof_y2, glob_dof_y3
    INTEGER              :: i,j
 
    
    m_mass      = 0.
    node_vol    = 0. !< @todo: this is not nice. 'calculate_node_vols' is already called by 'setup_stiffness_matrix'
                     !! we only need the parameters fi, gi and hi and the Area, however.
    !$omp parallel do private(A,m_mass_local,glob_dof_x1,glob_dof_x2,glob_dof_x3,glob_dof_y1,glob_dof_y2,glob_dof_y3) reduction(+:m_mass)
    DO i = 1, M_elems 
     
      CALL calculate_node_volumes(i, A)

      m_mass_local = 0.

      !--- fill the local mass matrix of the current element i ---
      m_mass_local(1,1) = 2. !N_1^1 
      m_mass_local(1,2) = 0.
      m_mass_local(1,3) = 1. !N_1*N_2
      m_mass_local(1,4) = 0.
      m_mass_local(1,5) = 1. !N_1*N_3
      m_mass_local(1,6) = 0.

      m_mass_local(2,1) = 0.
      m_mass_local(2,2) = 2. 
      m_mass_local(2,3) = 0.
      m_mass_local(2,4) = 1.
      m_mass_local(2,5) = 0.
      m_mass_local(2,6) = 1.

      m_mass_local(3,1) = 1.
      m_mass_local(3,2) = 0.
      m_mass_local(3,3) = 2. !N_2^2
      m_mass_local(3,4) = 0.
      m_mass_local(3,5) = 1. !N_2*N_3
      m_mass_local(3,6) = 0.

      m_mass_local(4,1) = 0.
      m_mass_local(4,2) = 1.
      m_mass_local(4,3) = 0.
      m_mass_local(4,4) = 2.
      m_mass_local(4,5) = 0.
      m_mass_local(4,6) = 1.

      m_mass_local(5,1) = 1.
      m_mass_local(5,2) = 0.
      m_mass_local(5,3) = 1.
      m_mass_local(5,4) = 0.
      m_mass_local(5,5) = 2. !N_3^2
      m_mass_local(5,6) = 0.

      m_mass_local(6,1) = 0.
      m_mass_local(6,2) = 1. 
      m_mass_local(6,3) = 0.
      m_mass_local(6,4) = 1.
      m_mass_local(6,5) = 0.
      m_mass_local(6,6) = 2.

      !--- multiply with remaining factors ---
      m_mass_local = (A/12.) * rho_solid * m_mass_local

      !--- calculate global mass matrix ---
      !--- global DOFs ---
      glob_dof_x1 = elems(i,1)*2 - 1
      glob_dof_y1 = elems(i,1)*2
      glob_dof_x2 = elems(i,2)*2 - 1
      glob_dof_y2 = elems(i,2)*2 
      glob_dof_x3 = elems(i,3)*2 - 1
      glob_dof_y3 = elems(i,3)*2


      !--- x-component of node 1 ---
      m_mass(glob_dof_x1, glob_dof_x1) = m_mass(glob_dof_x1, glob_dof_x1) + m_mass_local(1,1)
      
      m_mass(glob_dof_x1, glob_dof_y1) = m_mass(glob_dof_x1, glob_dof_y1) + m_mass_local(1,2)
                                           
      m_mass(glob_dof_x1, glob_dof_x2) = m_mass(glob_dof_x1, glob_dof_x2) + m_mass_local(1,3)
                                           
      m_mass(glob_dof_x1, glob_dof_y2) = m_mass(glob_dof_x1, glob_dof_y2) + m_mass_local(1,4)
                                     
      m_mass(glob_dof_x1, glob_dof_x3) = m_mass(glob_dof_x1, glob_dof_x3) + m_mass_local(1,5)
                                     
      m_mass(glob_dof_x1, glob_dof_y3) = m_mass(glob_dof_x1, glob_dof_y3) + m_mass_local(1,6)
                                    
                                       
      !--- y-component of node 1 ---
      m_mass(glob_dof_y1, glob_dof_x1) = m_mass(glob_dof_y1, glob_dof_x1) + m_mass_local(2,1)
                                          
      m_mass(glob_dof_y1, glob_dof_y1) = m_mass(glob_dof_y1, glob_dof_y1) + m_mass_local(2,2)
                                           
      m_mass(glob_dof_y1, glob_dof_x2) = m_mass(glob_dof_y1, glob_dof_x2) + m_mass_local(2,3)
                                           
      m_mass(glob_dof_y1, glob_dof_y2) = m_mass(glob_dof_y1, glob_dof_y2) + m_mass_local(2,4)
                                     
      m_mass(glob_dof_y1, glob_dof_x3) = m_mass(glob_dof_y1, glob_dof_x3) + m_mass_local(2,5)
                                     
      m_mass(glob_dof_y1, glob_dof_y3) = m_mass(glob_dof_y1, glob_dof_y3) + m_mass_local(2,6)
                                     

      !--- x-component of node 2 ---
      m_mass(glob_dof_x2, glob_dof_x1) = m_mass(glob_dof_x2, glob_dof_x1) + m_mass_local(3,1)
                                           
      m_mass(glob_dof_x2, glob_dof_y1) = m_mass(glob_dof_x2, glob_dof_y1) + m_mass_local(3,2)
                                           
      m_mass(glob_dof_x2, glob_dof_x2) = m_mass(glob_dof_x2, glob_dof_x2) + m_mass_local(3,3)
                                            
      m_mass(glob_dof_x2, glob_dof_y2) = m_mass(glob_dof_x2, glob_dof_y2) + m_mass_local(3,4)
                                    
      m_mass(glob_dof_x2, glob_dof_x3) = m_mass(glob_dof_x2, glob_dof_x3) + m_mass_local(3,5)
                                     
      m_mass(glob_dof_x2, glob_dof_y3) = m_mass(glob_dof_x2, glob_dof_y3) + m_mass_local(3,6)
                                     

      !--- y-component of node 2 ---
      m_mass(glob_dof_y2, glob_dof_x1) = m_mass(glob_dof_y2, glob_dof_x1) + m_mass_local(4,1)
                                           
      m_mass(glob_dof_y2, glob_dof_y1) = m_mass(glob_dof_y2, glob_dof_y1) + m_mass_local(4,2) 
                                           
      m_mass(glob_dof_y2, glob_dof_x2) = m_mass(glob_dof_y2, glob_dof_x2) + m_mass_local(4,3)
                                          
      m_mass(glob_dof_y2, glob_dof_y2) = m_mass(glob_dof_y2, glob_dof_y2) + m_mass_local(4,4)
                                     
      m_mass(glob_dof_y2, glob_dof_x3) = m_mass(glob_dof_y2, glob_dof_x3) + m_mass_local(4,5)
                                     
      m_mass(glob_dof_y2, glob_dof_y3) = m_mass(glob_dof_y2, glob_dof_y3) + m_mass_local(4,6)
                                     

      !--- x-component of node 3 ---
      m_mass(glob_dof_x3, glob_dof_x1) = m_mass(glob_dof_x3, glob_dof_x1) + m_mass_local(5,1)
                                           
      m_mass(glob_dof_x3, glob_dof_y1) = m_mass(glob_dof_x3, glob_dof_y1) + m_mass_local(5,2)
                                           
      m_mass(glob_dof_x3, glob_dof_x2) = m_mass(glob_dof_x3, glob_dof_x2) + m_mass_local(5,3)
                                           
      m_mass(glob_dof_x3, glob_dof_y2) = m_mass(glob_dof_x3, glob_dof_y2) + m_mass_local(5,4)                

      m_mass(glob_dof_x3, glob_dof_x3) = m_mass(glob_dof_x3, glob_dof_x3) + m_mass_local(5,5)
                                     
      m_mass(glob_dof_x3, glob_dof_y3) = m_mass(glob_dof_x3, glob_dof_y3) + m_mass_local(5,6)
                                     

      !--- y-component of node 3 ---
      m_mass(glob_dof_y3, glob_dof_x1) = m_mass(glob_dof_y3, glob_dof_x1) + m_mass_local(6,1)
                                           
      m_mass(glob_dof_y3, glob_dof_y1) = m_mass(glob_dof_y3, glob_dof_y1) + m_mass_local(6,2)
                                          
      m_mass(glob_dof_y3, glob_dof_x2) = m_mass(glob_dof_y3, glob_dof_x2) + m_mass_local(6,3)
                                           
      m_mass(glob_dof_y3, glob_dof_y2) = m_mass(glob_dof_y3, glob_dof_y2) + m_mass_local(6,4)
                                     
      m_mass(glob_dof_y3, glob_dof_x3) = m_mass(glob_dof_y3, glob_dof_x3) + m_mass_local(6,5)
                                     
      m_mass(glob_dof_y3, glob_dof_y3) = m_mass(glob_dof_y3, glob_dof_y3) + m_mass_local(6,6)

    END DO
    !$omp end parallel do


  END SUBROUTINE setup_mass_matrix




  SUBROUTINE setup_mass_matrix2

    IMPLICIT NONE

    REAL                 :: f1,f2,f3
    REAL                 :: g1,g2,g3
    REAL                 :: h1,h2,h3
    REAL                 :: x1, x2, x3, y1, y2, y3
    REAL                 :: A
    REAL, DIMENSION(6,6) :: m_mass_local

    INTEGER              :: glob_dof_x1, glob_dof_x2, glob_dof_x3, glob_dof_y1, glob_dof_y2, glob_dof_y3
    INTEGER              :: i,j
 
    
    m_mass      = 0.
    node_vol    = 0. !< @todo: this is not nice. 'calculate_node_vols' is already called by 'setup_stiffness_matrix'
                     !! we only need the parameters fi, gi and hi and the Area, however.
    !$omp parallel do private(A,f1,f2,f3,g1,g2,g3,h1,h2,h3,x1,x2,x3,y1,y2,y3,m_mass_local,glob_dof_x1,glob_dof_x2,glob_dof_x3,glob_dof_y1,glob_dof_y2,glob_dof_y3) reduction(+:m_mass)
    DO i = 1, M_elems
     
      CALL calculate_node_volumes(i, A, f1, f2, f3, g1, g2, g3, h1, h2, h3, x1, x2, x3, y1, y2, y3)

      m_mass_local = 0.

      !--- fill the local mass matrix of the current element i ---
      m_mass_local(1,1) = gq_int(g1,h1,f1,g1,h1,f1,x1,x2,x3,y1,y2,y3,A) !N_1^1 
      m_mass_local(1,2) = 0.
      m_mass_local(1,3) = gq_int(g1,h1,f1,g2,h2,f2,x1,x2,x3,y1,y2,y3,A) !N_1*N_2
      m_mass_local(1,4) = 0.
      m_mass_local(1,5) = gq_int(g1,h1,f1,g3,h3,f3,x1,x2,x3,y1,y2,y3,A) !N_1*N_3
      m_mass_local(1,6) = 0.

      m_mass_local(2,1) = 0.
      m_mass_local(2,2) = m_mass_local(1,1)
      m_mass_local(2,3) = 0.
      m_mass_local(2,4) = m_mass_local(1,3)
      m_mass_local(2,5) = 0.
      m_mass_local(2,6) = m_mass_local(1,5)

      m_mass_local(3,1) = m_mass_local(1,3)
      m_mass_local(3,2) = 0.
      m_mass_local(3,3) = gq_int(g2,h2,f2,g2,h2,f2,x1,x2,x3,y1,y2,y3,A) !N_2^2
      m_mass_local(3,4) = 0.
      m_mass_local(3,5) = gq_int(g1,h1,f1,g1,h1,f1,x1,x2,x3,y1,y2,y3,A) !N_2*N_3
      m_mass_local(3,6) = 0.

      m_mass_local(4,1) = 0.
      m_mass_local(4,2) = m_mass_local(3,1)
      m_mass_local(4,3) = 0.
      m_mass_local(4,4) = m_mass_local(3,3)
      m_mass_local(4,5) = 0.
      m_mass_local(4,6) = m_mass_local(3,5)

      m_mass_local(5,1) = m_mass_local(1,5)
      m_mass_local(5,2) = 0.
      m_mass_local(5,3) = m_mass_local(3,5)
      m_mass_local(5,4) = 0.
      m_mass_local(5,5) = gq_int(g1,h1,f1,g1,h1,f1,x1,x2,x3,y1,y2,y3,A) !N_3^2
      m_mass_local(5,6) = 0.

      m_mass_local(6,1) = 0.
      m_mass_local(6,2) = m_mass_local(5,1) 
      m_mass_local(6,3) = 0.
      m_mass_local(6,4) = m_mass_local(5,3)
      m_mass_local(6,5) = 0.
      m_mass_local(6,6) = m_mass_local(5,5)

      !--- multiply with remaining factors ---
      m_mass_local = 1/(4* A**2) * rho_solid * m_mass_local

      !--- calculate global mass matrix ---
      !--- global DOFs ---
      glob_dof_x1 = elems(i,1)*2 - 1
      glob_dof_y1 = elems(i,1)*2
      glob_dof_x2 = elems(i,2)*2 - 1
      glob_dof_y2 = elems(i,2)*2 
      glob_dof_x3 = elems(i,3)*2 - 1
      glob_dof_y3 = elems(i,3)*2


      !--- x-component of node 1 ---
      m_mass(glob_dof_x1, glob_dof_x1) = m_mass(glob_dof_x1, glob_dof_x1) + m_mass_local(1,1)
      
      m_mass(glob_dof_x1, glob_dof_y1) = m_mass(glob_dof_x1, glob_dof_y1) + m_mass_local(1,2)
                                           
      m_mass(glob_dof_x1, glob_dof_x2) = m_mass(glob_dof_x1, glob_dof_x2) + m_mass_local(1,3)
                                           
      m_mass(glob_dof_x1, glob_dof_y2) = m_mass(glob_dof_x1, glob_dof_y2) + m_mass_local(1,4)
                                     
      m_mass(glob_dof_x1, glob_dof_x3) = m_mass(glob_dof_x1, glob_dof_x3) + m_mass_local(1,5)
                                     
      m_mass(glob_dof_x1, glob_dof_y3) = m_mass(glob_dof_x1, glob_dof_y3) + m_mass_local(1,6)
                                    
                                       
      !--- y-component of node 1 ---
      m_mass(glob_dof_y1, glob_dof_x1) = m_mass(glob_dof_y1, glob_dof_x1) + m_mass_local(2,1)
                                          
      m_mass(glob_dof_y1, glob_dof_y1) = m_mass(glob_dof_y1, glob_dof_y1) + m_mass_local(2,2)
                                           
      m_mass(glob_dof_y1, glob_dof_x2) = m_mass(glob_dof_y1, glob_dof_x2) + m_mass_local(2,3)
                                           
      m_mass(glob_dof_y1, glob_dof_y2) = m_mass(glob_dof_y1, glob_dof_y2) + m_mass_local(2,4)
                                     
      m_mass(glob_dof_y1, glob_dof_x3) = m_mass(glob_dof_y1, glob_dof_x3) + m_mass_local(2,5)
                                     
      m_mass(glob_dof_y1, glob_dof_y3) = m_mass(glob_dof_y1, glob_dof_y3) + m_mass_local(2,6)
                                     

      !--- x-component of node 2 ---
      m_mass(glob_dof_x2, glob_dof_x1) = m_mass(glob_dof_x2, glob_dof_x1) + m_mass_local(3,1)
                                           
      m_mass(glob_dof_x2, glob_dof_y1) = m_mass(glob_dof_x2, glob_dof_y1) + m_mass_local(3,2)
                                           
      m_mass(glob_dof_x2, glob_dof_x2) = m_mass(glob_dof_x2, glob_dof_x2) + m_mass_local(3,3)
                                            
      m_mass(glob_dof_x2, glob_dof_y2) = m_mass(glob_dof_x2, glob_dof_y2) + m_mass_local(3,4)
                                    
      m_mass(glob_dof_x2, glob_dof_x3) = m_mass(glob_dof_x2, glob_dof_x3) + m_mass_local(3,5)
                                     
      m_mass(glob_dof_x2, glob_dof_y3) = m_mass(glob_dof_x2, glob_dof_y3) + m_mass_local(3,6)
                                     

      !--- y-component of node 2 ---
      m_mass(glob_dof_y2, glob_dof_x1) = m_mass(glob_dof_y2, glob_dof_x1) + m_mass_local(4,1)
                                           
      m_mass(glob_dof_y2, glob_dof_y1) = m_mass(glob_dof_y2, glob_dof_y1) + m_mass_local(4,2) 
                                           
      m_mass(glob_dof_y2, glob_dof_x2) = m_mass(glob_dof_y2, glob_dof_x2) + m_mass_local(4,3)
                                          
      m_mass(glob_dof_y2, glob_dof_y2) = m_mass(glob_dof_y2, glob_dof_y2) + m_mass_local(4,4)
                                     
      m_mass(glob_dof_y2, glob_dof_x3) = m_mass(glob_dof_y2, glob_dof_x3) + m_mass_local(4,5)
                                     
      m_mass(glob_dof_y2, glob_dof_y3) = m_mass(glob_dof_y2, glob_dof_y3) + m_mass_local(4,6)
                                     

      !--- x-component of node 3 ---
      m_mass(glob_dof_x3, glob_dof_x1) = m_mass(glob_dof_x3, glob_dof_x1) + m_mass_local(5,1)
                                           
      m_mass(glob_dof_x3, glob_dof_y1) = m_mass(glob_dof_x3, glob_dof_y1) + m_mass_local(5,2)
                                           
      m_mass(glob_dof_x3, glob_dof_x2) = m_mass(glob_dof_x3, glob_dof_x2) + m_mass_local(5,3)
                                           
      m_mass(glob_dof_x3, glob_dof_y2) = m_mass(glob_dof_x3, glob_dof_y2) + m_mass_local(5,4)                

      m_mass(glob_dof_x3, glob_dof_x3) = m_mass(glob_dof_x3, glob_dof_x3) + m_mass_local(5,5)
                                     
      m_mass(glob_dof_x3, glob_dof_y3) = m_mass(glob_dof_x3, glob_dof_y3) + m_mass_local(5,6)
                                     

      !--- y-component of node 3 ---
      m_mass(glob_dof_y3, glob_dof_x1) = m_mass(glob_dof_y3, glob_dof_x1) + m_mass_local(6,1)
                                           
      m_mass(glob_dof_y3, glob_dof_y1) = m_mass(glob_dof_y3, glob_dof_y1) + m_mass_local(6,2)
                                          
      m_mass(glob_dof_y3, glob_dof_x2) = m_mass(glob_dof_y3, glob_dof_x2) + m_mass_local(6,3)
                                           
      m_mass(glob_dof_y3, glob_dof_y2) = m_mass(glob_dof_y3, glob_dof_y2) + m_mass_local(6,4)
                                     
      m_mass(glob_dof_y3, glob_dof_x3) = m_mass(glob_dof_y3, glob_dof_x3) + m_mass_local(6,5)
                                     
      m_mass(glob_dof_y3, glob_dof_y3) = m_mass(glob_dof_y3, glob_dof_y3) + m_mass_local(6,6)

    END DO
    !$omp end parallel do


  END SUBROUTINE setup_mass_matrix2



  !> function to integrate the function 'integrand' over quadrature points of the specified order N
  FUNCTION gq_int(b1, c1, d1, b2, c2, d2, x1, x2, x3, y1, y2, y3, A) RESULT(fn_val)

    IMPLICIT NONE

    REAL   , INTENT(IN) :: b1, c1, d1, b2, c2, d2
    REAL   , INTENT(IN) :: x1, x2, x3, y1, y2, y3, A
    REAL                :: fn_val
    INTEGER, PARAMETER  :: N    = 3 !< order of Gaussian quadrature (only assigned during first call in FORTRAN)
    INTEGER, PARAMETER  :: dims = 4 !< row dimension of wqp array that store qp's (hardcoded for now) look into 'qp_lookup'
                                    !! if wanting to change
    REAL   , DIMENSION(dims,3) :: wqp
    INTEGER             :: i,j
    REAL                :: x,y


    !--- first retrieve the quadrature points' coordinate in the reference triangle and the corresponding weights ---
    CALL qp_lookup(N, dims, wqp)

    !--- loop over all Gaussian quadrature points to integrate ---
    fn_val = 0.
    !$omp parallel do private(x,y) reduction(+:fn_val)
    DO i = 1, dims 
      x = x1*(1 - wqp(i,1) -  wqp(i,2)) + x2*wqp(i,1) + x3*wqp(i,2)
      y = y1*(1 - wqp(i,1) -  wqp(i,2)) + y2*wqp(i,1) + y3*wqp(i,2)
      fn_val = fn_val + integrand(b1, c1, d1, b2, c2, d2,x, y)*wqp(i,3)
    END DO
    !$omp end parallel do
    fn_val = A*fn_val

    RETURN

  END FUNCTION gq_int

  !> function that returns the integrand of one matrix entry for the mass_matrix
  FUNCTION integrand(b1, c1, d1, b2, c2, d2, x, y) RESULT(fn_val)

  IMPLICIT NONE

  REAL, INTENT(IN) :: b1, c1, d1, b2, c2, d2
  REAL, INTENT(IN) :: x, y
  REAL             :: fn_val

  fn_val = (x*b1 + y*c1 + d1)*(x*b2 + y*c2 + d2)

  RETURN


  END FUNCTION


  SUBROUTINE qp_lookup(N, dims, wqp)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, dims
  REAL   , INTENT(OUT):: wqp(dims,3)

  integer          :: i

  IF (N .EQ. 1) THEN
    wqp = TRANSPOSE(RESHAPE( (/ 0.33333333333333, 0.33333333333333, 1.00000000000000 /), (/SIZE(wqp,2),SIZE(wqp,1)/) ) )

  ELSEIF (N .EQ. 2) THEN
    wqp = TRANSPOSE(RESHAPE( (/ 0.16666666666667, 0.16666666666667, 0.33333333333333, &
                                0.16666666666667, 0.66666666666667, 0.33333333333333, &
                                0.66666666666667, 0.16666666666667, 0.33333333333333 /), (/SIZE(wqp,2),SIZE(wqp,1)/) ) )
  ELSEIF (N .EQ. 3) THEN
    wqp = TRANSPOSE(RESHAPE( (/ 0.33333333333333, 0.33333333333333,-0.56250000000000, &
                                0.20000000000000, 0.20000000000000, 0.52083333333333, &
                                0.20000000000000, 0.60000000000000, 0.52083333333333, &
                                0.60000000000000, 0.20000000000000, 0.52083333333333 /), (/SIZE(wqp,2),SIZE(wqp,1)/) ) )
  ELSEIF (N .EQ. 4) THEN
    wqp = TRANSPOSE(RESHAPE( (/ 0.44594849091597, 0.44594849091597, 0.22338158967801, &
                                0.44594849091597, 0.10810301816807, 0.22338158967801, &
                                0.10810301816807, 0.44594849091597, 0.22338158967801, &
                                0.09157621350977, 0.09157621350977, 0.10995174365532, &
                                0.09157621350977, 0.81684757298046, 0.10995174365532, &
                                0.81684757298046, 0.09157621350977, 0.10995174365532 /), (/SIZE(wqp,2),SIZE(wqp,1)/) ) )
  ELSEIF (N .EQ. 5) THEN
    wqp = TRANSPOSE(RESHAPE( (/ 0.33333333333333, 0.33333333333333, 0.22500000000000, &
                                0.47014206410511, 0.47014206410511, 0.13239415278851, &
                                0.47014206410511, 0.05971587178977, 0.13239415278851, &
                                0.05971587178977, 0.47014206410511, 0.13239415278851, &
                                0.10128650732346, 0.10128650732346, 0.12593918054483, &
                                0.10128650732346, 0.79742698535309, 0.12593918054483, &
                                0.79742698535309, 0.10128650732346, 0.12593918054483 /), (/SIZE(wqp,2),SIZE(wqp,1)/) ) )
  ELSEIF (N .EQ. 6) THEN
    wqp = TRANSPOSE(RESHAPE( (/ 0.24928674517091, 0.24928674517091, 0.11678627572638, &
                                0.24928674517091, 0.50142650965818, 0.11678627572638, &
                                0.50142650965818, 0.24928674517091, 0.11678627572638, &
                                0.06308901449150, 0.06308901449150, 0.05084490637021, &
                                0.06308901449150, 0.87382197101700, 0.05084490637021, &
                                0.87382197101700, 0.06308901449150, 0.05084490637021, &
                                0.31035245103378, 0.63650249912140, 0.08285107561837, &
                                0.63650249912140, 0.05314504984482, 0.08285107561837, &
                                0.05314504984482, 0.31035245103378, 0.08285107561837, &
                                0.63650249912140, 0.31035245103378, 0.08285107561837, &
                                0.31035245103378, 0.05314504984482, 0.08285107561837, &
                                0.05314504984482, 0.63650249912140, 0.08285107561837 /), (/SIZE(wqp,2),SIZE(wqp,1)/) ) )
   ELSE
     write(*,*) 'ERROR: Bad input N for quadrature order!'
     CALL MPI_FINALIZE(merror)
     STOP
   END IF

  END SUBROUTINE






  !> subroutine to compute the Lagrangian force (this will eventually be replaced by an FEM solver)
  SUBROUTINE compute_force
    
    IMPLICIT NONE

    INTEGER            :: i,j
    INTEGER            :: i_odd,i_even,j_slow,j_fast
    
    fb = 0.

    IF (fem_yes .AND. dimens .EQ. 2) THEN
    !$omp parallel do private(j,i_odd,i_even,j_slow,j_fast) reduction(+:fb)
    FORCE_LOOP: DO i = 1, M_bound ! loops over force components
      j_slow = 0
      IF ( ANY( i == NINT(ebcs(:,1) ) ) ) cycle FORCE_LOOP
      DO j = 1, 2*M_bound   ! loops over displacement components
      !--- compute the force according to Cauchy's stress equation --- 
                         i_odd   = 2*i - 1
                         i_even  = 2*i
                         j_slow  = MOD( j   ,2)+j_slow
                         j_fast  = MOD((j+1),2)+1

                         fb(i,1) = fb(i,1) + k_stiff(i_odd ,j)*db( j_slow ,j_fast ) &
                                           + m_mass(i_odd  ,j)*(ub(j_slow, j_fast) - container(j_slow,j_fast))/(dtime*(aRK(substep)+bRK(substep)))
                                           !+ c_damp(i_odd  ,j)*ub( j_slow, j_fast ) 
                         fb(i,2) = fb(i,2) + k_stiff(i_even,j)*db( j_slow ,j_fast ) &
                                           + m_mass(i_even ,j)*(ub(j_slow, j_fast) - container(j_slow,j_fast))/(dtime*(aRK(substep)+bRK(substep)))
                                           !+ c_damp(i_even ,j)*ub( j_slow, j_fast ) 
        IF (dimens == 3) write(*,*)'    WARNING: 3D FEM force not supported' 
      END DO
    END DO FORCE_LOOP
    !$omp end parallel do
    ELSE IF (.NOT. fem_yes .AND. dimens .EQ. 3) THEN
    
    SPRING_FORCE_LOOP: DO i = 1, M_bound
      IF ( ANY( i == NINT(ebcs(:,1) ) ) ) cycle SPRING_FORCE_LOOP
      IF (i .LE. M_ebcs) THEN
        fb(i,1) = 1000.*k_x * db(i,1)
        fb(i,2) = 1000.*k_y * db(i,2)
        fb(i,3) =       k_z * db(i,3)
      END IF
      !--- oscillator model
      fb(i,1) = k_x * db(i,1)
      fb(i,2) = k_y * db(i,2)
      fb(i,3) = k_z * db(i,3)

    END DO SPRING_FORCE_LOOP

    END IF

    !write(*,*) 'in process',rank,'fb=',fb
 
  END SUBROUTINE compute_force




  !> function to compute a one-dimensional discrete delta function
  !! @note can the return value be of other type than a scalar? (array, tensor). If not, consider a subroutine instead.
  FUNCTION discrete_delta(distance) RESULT(delta)
   
    IMPLICIT NONE
    
    REAL, INTENT(IN)     :: distance
    REAL                 :: delta, pi

    
    pi = 2.*ABS(ACOS(0.))

    IF      (ddf_type == 1) THEN
      IF      (ABS(distance) .LE. reach) THEN
        delta = ( 1. + COS( pi*ABS(distance)/reach ) )/( 2.*reach )
      ELSE IF (distance == L1*L2*L3) THEN
        delta = 0.
      ELSE
        delta = 0.
        WRITE(*,*)'     WARNING: zero kernel entry created. Distance is past reach.'
      END IF
    ELSE IF (ddf_type == 2) THEN
      IF      (ABS(distance) .LE. reach) THEN
        delta = ( reach - ABS(distance) )/( (reach**2) )
      ELSE IF (distance == L1*L2*L3) THEN
        delta = 0.
      ELSE
        delta = 0.
        WRITE(*,*)'     WARNING: zero kernel entry created. Distance is past reach.'
      END IF
    ELSE IF (ddf_type == 3) THEN
      IF      (ABS(distance) .LE. (reach/2.)) THEN
        delta = ( 1. - ( distance/(reach/2.) )**2 )/(reach/2.)
      ELSE IF ((ABS(distance) .GT. (reach/2.)) .AND. (ABS(distance) .LE. reach)) THEN
        delta = ( 2. - ABS(3.*distance/(reach/2.)) + (distance/(reach/2.))**2.)/(reach/2.)
      ELSE IF (distance == L1*L2*L3) THEN
        delta = 0.
      ELSE
        delta = 0.
        WRITE(*,*)'     WARNING: zero kernel entry created. Distance is past reach.'
      END IF
    ELSE IF (ddf_type == 4) THEN
      IF      (ABS(distance) .LE. (reach/2.)) THEN
        delta = ( 3. - 2.*ABS(distance)/(reach/2.) + SQRT( 1. + 4.*ABS(distance)/(reach/2.) - 4.*(ABS(distance)/(reach/2.))**2.))/( 4.*reach)
      ELSE IF ((ABS(distance) .GT. (reach/2.)).AND. (ABS(distance) .LE. reach)) THEN
        delta = ( 5. - 2.*ABS(distance)/(reach/2.) - SQRT(-7. +12.*ABS(distance)/(reach/2.) - 4.*(ABS(distance)/(reach/2.))**2.))/( 4.*reach)
      ELSE IF (distance == L1*L2*L3) THEN
        delta = 0.
      ELSE
        delta = 0.
        WRITE(*,*)'     WARNING: zero kernel entry created. Distance is past reach.'
      END IF
    END IF
       


  END FUNCTION discrete_delta


END MODULE mod_ibm
