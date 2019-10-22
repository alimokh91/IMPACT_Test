!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!*************************************************************************************************************

!> module containing the subroutine timeintegration, which is basically the central part of the DNS.
!! It uses modules mod_dims, mod_vars, mod_exchange, mod_diff, mod_inout, mod_coeffs, mod_lib, mod_test, mod_particles,
!! mod_solvers, mod_rhs
MODULE mod_timeint
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff
  USE mod_inout
  USE mod_coeffs
  USE mod_lib
  USE mod_test
  USE mod_solvers
  USE mod_rhs
  USE MPI

  PRIVATE
  
  PUBLIC timeintegration
  
  CONTAINS
  
!!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> public subroutine that forms the center of the DNS simulation.
  SUBROUTINE timeintegration
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  REAL                   ::  min_dx, Re_Ko 
 
  
  !--- Startvektor -------------------------------------------------------------------------------------------
  ! Hier bereits notwendig, um weitere Konzentrationen zuschalten zu koennen (BC werden aber weiter unten initialisiert!).
                         CALL initial_conditions_vel
  
  !--- diverse Files <F6>ffnen ----------------------------------------------------------------------------------
  IF (dtime_out_scal /= 0.) CALL open_stats
  IF (dtime_out_kalm /= 0.) CALL open_kalman
  
  IF (restart == 0) THEN
     time          = time_start
     time_out_vect = time_start
     time_out_scal = time_start
     time_out_kalm = time_start
     
     dtime         = 0.
     timestep      = 0
     
     new_dtime      = .TRUE.
     write_out_vect = .TRUE.
     write_out_scal = .TRUE.
     write_out_kalm = .TRUE.
     
     write_count = 0
     
     IF (dtime_out_vect == 0.) write_out_vect = .FALSE.
     IF (dtime_out_scal == 0.) write_out_scal = .FALSE.
     IF (dtime_out_kalm == 0.) write_out_kalm = .FALSE.
  ELSE
                               CALL read_restart
     IF (dtime_out_scal /= 0.) CALL read_restart_stats
     IF (dtime_out_kalm /= 0.) CALL read_restart_kalman
     
     IF (rank == 0 .AND. write_stout_yes) THEN
        WRITE(*,'(a,1E13.5)') '             time =', time
        WRITE(*,'(a,1i5   )') '         timestep =', timestep
        
        WRITE(*,'(a,1E13.5)') '            dtime =', dtime
        WRITE(*,'(a,1E13.5)') '    time_out_vect =', time_out_vect
        WRITE(*,'(a,1E13.5)') '    time_out_scal =', time_out_scal
        WRITE(*,'(a,1E13.5)') '    time_out_kalm =', time_out_kalm
        
        WRITE(*,'(a,1L2   )') '        new_dtime =', new_dtime
        WRITE(*,'(a,1L2   )') '   write_out_vect =', write_out_vect
        WRITE(*,'(a,1L2   )') '   write_out_scal =', write_out_scal
        WRITE(*,'(a,1L2   )') '   write_out_kalm =', write_out_kalm
        
        WRITE(*,'(a,1i8   )') '      write_count =', write_count
     END IF
  END IF
  
  timestep_old  = timestep
  dtime_old     = dtime
  dtime_average = 0.
  finish_yes    = .FALSE.
  
  !--- Null-Raeume bestimmen ---------------------------------------------------------------------------------
  ! Steht hier, weil Korrekturvektor "th" nach "configuration" erst alloziert und bestimmt werden muss
  ! ("initial_conditions_" werden danach als nächstes gerufen, s.o.)
  ! Alternative: eigene Subroutine für "th" kreieren ...
  IF (nullspace_yes) THEN
     CALL get_stencil_transp
     CALL get_nullspace
     CALL get_stencil ! TEST!!! Unschoen! besser zu impact.f90 packen und Reihenfolge abaendern ...
  END IF
  
  !--- RB initialisieren -------------------------------------------------------------------------------------
  CALL init_BC
  
  !--- Divergenz-Freiheit testen -----------------------------------------------------------------------------
  CALL test_divergence

  !--- File fuer zwischenzeitlichen Abbruch der Zeitintegration neu erstellen --------------------------------
  IF (rank == 0) THEN
     OPEN (10,FILE='send_signal.txt',STATUS='UNKNOWN')
     WRITE(10,'(a)') '0'
     WRITE(10,*) n_timesteps
     WRITE(10,*) time_end
     CALL flush(10)
     CLOSE(10)
  END IF
 
  !--- Compute the "Kolmogorov-Re" number (based on the smallest scale) ---------------------------------------
  min_dx = MIN(MINVAL(dy1p), MINVAL(dy2p), MINVAL(dy3p))
  Re_Ko  = rho_fluid*U_ref*min_dx/mu_fluid
 
  
  IF (rank == 0 .AND. write_stout_yes) THEN
     WRITE(*,'(a)')
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)')
     WRITE(*,'(a)') '---------------------------- START TIME-INTEGRATION -----------------------------'
     WRITE(*,'(a)')
     WRITE(*,'(a,E12.5)')                 '                     Re =',Re
     WRITE(*,'(a,E12.5)')                 '                     Re_Ko =',Re_Ko
     WRITE(*,'(a,i4,a,i4,a,i4)')          '          box resolution:',M1,' x',M2,' x',M3
     WRITE(*,'(a,E12.5,a,E12.5,a,E12.5)') '          box dimension :',L1,' x',L2,' x',L3
     WRITE(*,'(a,E12.5)')                 '                   epsU =',epsU
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)')
  END IF
  
  
  !--- Zeitmessung starten -----------------------------------------------------------------------------------
  IF (rank == 0) THEN
     CALL DATE_AND_TIME(values=ctime)
     day  = ctime(3)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     elatime = msec+1000*(sec+60*(minu+60*hour))
     OPEN(99,FILE='test_wallclocktime_restart'//restart_char//'.txt',STATUS='UNKNOWN')
     WRITE(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Begin time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &      ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     CALL flush(99)
  END IF
  
  
  !--- Ausschreiben ------------------------------------------------------------------------------------------
  IF (write_xdmf_yes .AND. write_out_vect) CALL write_xdmf_xml ! bbecsek
  IF (write_out_scal) CALL compute_stats
  IF (write_out_kalm .and. time.eq.time_start) CALL compute_kalman
  IF (write_out_vect) CALL write_fields
  !===========================================================================================================
  
  !--- ghost cell update -------------------------------------------------------------------------------
  CALL exchange_all_all(.TRUE.,vel)

  !--- interpolate advection velocity + update ghost cells ---------------------------------------------
  ! vel(:,:,:,i) --> worki(:,:,:)
  CALL interpolate_vel(.FALSE.) ! TEST!!! Wurde teilweise schon bei Zeitschritt-Bestimmung erledigt!

  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  timeint: DO
     
     CALL get_dtime
     
     IF (rank == 0 .AND. write_stout_yes) THEN
        IF (timeint_mode == 0) THEN
           WRITE(*,'(a)')
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a)') '================================================================================='
        END IF
        WRITE(*,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,' ; time =',time,' ; dtime =',dtime
        IF (timeint_mode == 0) THEN
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a)') '================================================================================='
        END IF
     END IF
     
     IF (rank == 0 .AND. log_iteration_yes) THEN
        OPEN(10,FILE='log_iterations.txt', STATUS='UNKNOWN')
        WRITE(10,'(a)')
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; dtime =',dtime
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a)') '================================================================================='
     END IF
     
     !========================================================================================================
     

     DO substep = 1, RK_steps
  
        
        IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) THEN
           WRITE(*,'(a)') 
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           WRITE(*,'(a)') '================================================================================='
        END IF
        IF (rank == 0 .AND. log_iteration_yes) THEN
           WRITE(10,'(a)') 
           WRITE(10,'(a)') '================================================================================='
           WRITE(10,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           WRITE(10,'(a)') '================================================================================='
        END IF
        
        !--- Zeit --------------------------------------------------------------------------------------------
        IF (substep == 1) subtime = time + dtime* aRK(1)
        IF (substep == 2) subtime = time + dtime*(aRK(1)+aRK(2)+bRK(2))
        IF (substep == 3) subtime = time + dtime
        
        !--- rhs (ggf. Neumann-RB überschreiben) -------------------------------------------------------------
        ! Muss vor Konzentrationen kommen, weil
        !  - bcii für die no-flux-RB verwendet werden,
        !  - die Konzentrationen die Eddy-Viscosity des Geschwindigkeitsfeldes benoetigen.
        CALL rhs_vel
        
        !--- Helmholtz-Multiplikator -------------------------------------------------------------------------
        multL = thetaL*(aRK(substep)+bRK(substep))*dtime / Re
        
        !--- Umskalieren (Effizienz, initial guess) ----------------------------------------------------------
        IF (.NOT. init_pre(substep)) pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) * (aRK(substep)+bRK(substep)) * dtime
        
        !--- Löser -------------------------------------------------------------------------------------------
        IF (timeint_mode == 1 .OR. thetaL == 1.) THEN
           CALL explicit
        ELSE
           IF (twostep_yes) THEN
              CALL twostep
           ELSE
              CALL outer_iteration
           END IF
        END IF
        
        !--- physikalischer Druck ----------------------------------------------------------------------------
        pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) / (aRK(substep)+bRK(substep)) / dtime
        
        !--- Undefinierte Ecken / Kanten auffüllen -----------------------------------------------------------
        CALL fill_corners(pre)

        !--- ghost cell update (fuer RHS) --------------------------------------------------------------------
        CALL exchange_all_all(.TRUE.,vel)

        !--- interpolate advection velocity + update ghost cells ---------------------------------------------
        ! vel(:,:,:,i) --> worki(:,:,:)
        CALL interpolate_vel(.FALSE.) ! TEST!!! Wurde teilweise schon bei Zeitschritt-Bestimmung erledigt!

     END DO
     
     !========================================================================================================
     timestep = timestep + 1
     time     = time + dtime
     
     !--- send_signal.txt lesen ------------------------------------------------------------------------------
     CALL check_signal
     
     !--- Druck-Niveau festhalten ----------------------------------------------------------------------------
     CALL level_pressure
     
     !--- Ausschreiben ---------------------------------------------------------------------------------------
     IF (write_xdmf_yes .AND. write_out_vect) CALL write_xdmf_xml ! bbecsek
     IF (write_out_scal) CALL compute_stats
     IF (write_out_kalm) CALL compute_kalman
     IF (write_out_vect) CALL write_fields
     
     !--------------------------------------------------------------------------------------------------------
     IF (rank == 0 .AND. log_iteration_yes) CLOSE(10)
     
     IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) WRITE(*,*)
     IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) WRITE(*,*)
     
     CALL MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror) ! notwendig fuer "check_alarm"
     IF (finish_yes) EXIT
     
  END DO timeint
  !===========================================================================================================
  
  !--- Zeitmessung beenden -----------------------------------------------------------------------------------
  IF (rank == 0) THEN
     CALL DATE_AND_TIME(values=ctime)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     IF (ctime(3) /= day) THEN
        ! Anmerkung: Gilt nur für Jobs <= 24h
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime + 24*60*60*1000
     ELSE
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime
     END IF
     
     WRITE(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Finish time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &       ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     WRITE(99,'(a,E13.5)') 'elapsed time [sec]', REAL(elatime)/1000.
     CLOSE(99)
  END IF
  
  !--- Restart schreiben -------------------------------------------------------------------------------------
  restart = restart + 1
  CALL write_restart
  IF (dtime_out_scal /= 0.) CALL write_restart_stats
  IF (dtime_out_kalm /= 0.) CALL write_restart_kalman
  
  !--- Iterationsstatistiken auswerten -----------------------------------------------------------------------
  CALL iteration_stats
  
  !--- diverse Files schliessen ------------------------------------------------------------------------------
  IF (dtime_out_scal /= 0.) CALL close_stats
  IF (dtime_out_kalm /= 0.) CALL close_kalman
 
  !--- link all XDMF files together --------------------------------------------------------------------------
  IF (write_xdmf_yes) CALL write_xdmf_timecollection ! bbecsek

  END SUBROUTINE timeintegration
  
  
  
END MODULE mod_timeint
