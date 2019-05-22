!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!* modified by Barna Becsek, ARTORG CVE, University of Bern (barna.becsek@artorg.unibe.ch)                   *
!* October 2014                                                                                              *
!*************************************************************************************************************


!> module holding subroutines for value exchange between blocks (ghost cells). It uses modules mod_dims, mod_vars
!! @todo: there is two practically identical subroutines exchange and exchange2, with the latter containing variable
!!        interval limits. Could one be obsolete? If so, perhaps every call to the obsolete routine could be updated
!!        in the entire code and the obsolete routine could be removed.
MODULE mod_exchange
  
  
  USE mod_dims
  USE mod_vars
  !USE mpi

  PRIVATE
  
  PUBLIC exchange, exchange2, exchange_all_all, exchange_relax, exchange_part
  
  INCLUDE 'mpif.h'
  
  INTEGER                ::  status(MPI_STATUS_SIZE) ! TEST!!! Warum steht das eigentlich hier lokal? gleiches gilt auch fuer die anderen Module ...
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  !> subroutine for exchanging data between ghost cells.
  !! @note: this routine might be deprecated since exchange2 is identical with added adjustment of interval limits.
  SUBROUTINE exchange(dir,vel_dir,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in)    ::  dir      !< spacial direction
  INTEGER, INTENT(in)    ::  vel_dir  !< direction of velocity component
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) !< field that is exchanged 
  
  REAL                   ::  ghost1LR( 1:b1U,1:N2,1:N3) ! TEST!!! Feld global anlegen und all_all-Variante loeschen? Performance ...
  REAL                   ::  ghost1LS( 1:b1U,1:N2,1:N3)
  REAL                   ::  ghost1UR(b1L:-1,1:N2,1:N3)
  REAL                   ::  ghost1US(b1L:-1,1:N2,1:N3)
  
  REAL                   ::  ghost2LR(1:N1, 1:b2U,1:N3)
  REAL                   ::  ghost2LS(1:N1, 1:b2U,1:N3)
  REAL                   ::  ghost2UR(1:N1,b2L:-1,1:N3)
  REAL                   ::  ghost2US(1:N1,b2L:-1,1:N3)
  
  REAL                   ::  ghost3LR(1:N1,1:N2, 1:b3U)
  REAL                   ::  ghost3LS(1:N1,1:N2, 1:b3U)
  REAL                   ::  ghost3UR(1:N1,1:N2,b3L:-1)
  REAL                   ::  ghost3US(1:N1,1:N2,b3L:-1)
  
  INTEGER                ::  length1L, length1U
  INTEGER                ::  length2L, length2U
  INTEGER                ::  length3L, length3U
  
  INTEGER                ::  i, j, k
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  IF (dir == 1) THEN
     
     length1L = N2*N3*ABS(b1L)
     length1U = N2*N3*ABS(b1U)
     
     IF (BC_1L == 0) CALL MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
     IF (BC_1U == 0) CALL MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
     
     IF (BC_1U == 0) THEN
        DO i = b1L, -1
           ghost1US(i,1:N2,1:N3) = phi((N1p+1+i),1:N2,1:N3)
        END DO
     END IF
     IF (BC_1L == 0) THEN
        DO i = 1, b1U
           ghost1LS(i,1:N2,1:N3) = phi((S1p-1+i),1:N2,1:N3)
        END DO
     END IF
     
     IF (BC_1U == 0) CALL MPI_SEND(ghost1US,length1L,MPI_REAL8,rank1U,1,COMM_CART,merror)
     IF (BC_1L == 0) CALL MPI_SEND(ghost1LS,length1U,MPI_REAL8,rank1L,2,COMM_CART,merror)
     
     IF (BC_1L == 0) CALL MPI_WAIT(req1L,status,merror)
     IF (BC_1U == 0) CALL MPI_WAIT(req1U,status,merror)
     
     IF (BC_1L == 0) CALL pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     IF (BC_1U == 0) CALL pseudocall(ghost1LR)
     
     IF (BC_1L == 0) THEN
        DO i = b1L, -1
           phi((S1p+i),1:N2,1:N3) = ghost1UR(i,1:N2,1:N3)
        END DO
     END IF
     IF (BC_1U == 0) THEN
        DO i = 1, b1U
           phi((N1p+i),1:N2,1:N3) = ghost1LR(i,1:N2,1:N3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (BC_1L == -1) THEN
        DO i = b1L, -1
           phi((S1p+i),1:N2,1:N3) = phi((N1p+1+i),1:N2,1:N3)
        END DO
        DO i = 1, b1U
           phi((N1p+i),1:N2,1:N3) = phi((S1p-1+i),1:N2,1:N3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (vel_dir == dir) THEN
        IF (BC_1L .GT. 0) phi( b1L  : -1     ,1:N2,1:N3) = 0.
        IF (BC_1U .GT. 0) phi((N1+1):(N1+b1U),1:N2,1:N3) = 0.
        
        IF (BC_1L  == -2) phi( b1L  :  0     ,1:N2,1:N3) = 0.
        IF (BC_1U  == -2) phi( N1   :(N1+b1U),1:N2,1:N3) = 0.
     ELSE
        IF (BC_1L .GT. 0 .OR. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3) = 0.
        IF (BC_1U .GT. 0 .OR. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3) = 0.
     END IF
     
  END IF
  
  !======================================================================================================================
  
  IF (dir == 2) THEN
     
     length2L = N1*N3*ABS(b2L)
     length2U = N1*N3*ABS(b2U)
     
     IF (BC_2L == 0) CALL MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
     IF (BC_2U == 0) CALL MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
     
     IF (BC_2U == 0) THEN
        DO j = b2L, -1
           ghost2US(1:N1,j,1:N3) = phi(1:N1,(N2p+1+j),1:N3)
        END DO
     END IF
     IF (BC_2L == 0) THEN
        DO j = 1, b2U
           ghost2LS(1:N1,j,1:N3) = phi(1:N1,(S2p-1+j),1:N3)
        END DO
     END IF
     
     IF (BC_2U == 0) CALL MPI_SEND(ghost2US,length2L,MPI_REAL8,rank2U,3,COMM_CART,merror)
     IF (BC_2L == 0) CALL MPI_SEND(ghost2LS,length2U,MPI_REAL8,rank2L,4,COMM_CART,merror)
     
     IF (BC_2L == 0) CALL MPI_WAIT(req2L,status,merror)
     IF (BC_2U == 0) CALL MPI_WAIT(req2U,status,merror)
     
     IF (BC_2L == 0) CALL pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     IF (BC_2U == 0) CALL pseudocall(ghost2LR)
     
     IF (BC_2L == 0) THEN
        DO j = b2L, -1
           phi(1:N1,(S2p+j),1:N3) = ghost2UR(1:N1,j,1:N3)
        END DO
     END IF
     IF (BC_2U == 0) THEN
        DO j = 1, b2U
           phi(1:N1,(N2p+j),1:N3) = ghost2LR(1:N1,j,1:N3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (BC_2L == -1) THEN
        DO j = b2L, -1
           phi(1:N1,(S2p+j),1:N3) = phi(1:N1,(N2p+1+j),1:N3)
        END DO
        DO j = 1, b2U
           phi(1:N1,(N2p+j),1:N3) = phi(1:N1,(S2p-1+j),1:N3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (vel_dir == dir) THEN
        IF (BC_2L .GT. 0) phi(1:N1, b2L  : -1     ,1:N3) = 0.
        IF (BC_2U .GT. 0) phi(1:N1,(N2+1):(N2+b2U),1:N3) = 0.
        
        IF (BC_2L  == -2) phi(1:N1, b2L  :  0     ,1:N3) = 0.
        IF (BC_2U  == -2) phi(1:N1, N2   :(N2+b2U),1:N3) = 0.
     ELSE
        IF (BC_2L .GT. 0 .OR. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3) = 0.
        IF (BC_2U .GT. 0 .OR. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3) = 0.
     END IF
     
  END IF
  
  !======================================================================================================================
  
  IF (dir == 3 .AND. dimens == 3) THEN
     
     length3L = N1*N2*ABS(b3L)
     length3U = N1*N2*ABS(b3U)
     
     IF (BC_3L == 0) CALL MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
     IF (BC_3U == 0) CALL MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
     
     IF (BC_3U == 0) THEN
        DO k = b3L, -1
           ghost3US(1:N1,1:N2,k) = phi(1:N1,1:N2,(N3p+1+k))
        END DO
     END IF
     IF (BC_3L == 0) THEN
        DO k = 1, b3U
           ghost3LS(1:N1,1:N2,k) = phi(1:N1,1:N2,(S3p-1+k))
        END DO
     END IF
     
     IF (BC_3U == 0) CALL MPI_SEND(ghost3US,length3L,MPI_REAL8,rank3U,5,COMM_CART,merror)
     IF (BC_3L == 0) CALL MPI_SEND(ghost3LS,length3U,MPI_REAL8,rank3L,6,COMM_CART,merror)
     
     IF (BC_3L == 0) CALL MPI_WAIT(req3L,status,merror)
     IF (BC_3U == 0) CALL MPI_WAIT(req3U,status,merror)
     
     IF (BC_3L == 0) CALL pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     IF (BC_3U == 0) CALL pseudocall(ghost3LR)
     
     IF (BC_3L == 0) THEN
        DO k = b3L, -1
           phi(1:N1,1:N2,(S3p+k)) = ghost3UR(1:N1,1:N2,k)
        END DO
     END IF
     IF (BC_3U == 0) THEN
        DO k = 1, b3U
           phi(1:N1,1:N2,(N3p+k)) = ghost3LR(1:N1,1:N2,k)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (BC_3L == -1) THEN
        DO k = b3L, -1
           phi(1:N1,1:N2,(S3p+k)) = phi(1:N1,1:N2,(N3p+1+k))
        END DO
        DO k = 1, b3U
           phi(1:N1,1:N2,(N3p+k)) = phi(1:N1,1:N2,(S3p-1+k))
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (vel_dir == dir) THEN
        IF (BC_3L .GT. 0) phi(1:N1,1:N2, b3L  : -1     ) = 0.
        IF (BC_3U .GT. 0) phi(1:N1,1:N2,(N3+1):(N3+b3U)) = 0.
        
        IF (BC_3L  == -2) phi(1:N1,1:N2, b3L  :  0     ) = 0.
        IF (BC_3U  == -2) phi(1:N1,1:N2, N3   :(N3+b3U)) = 0.
     ELSE
        IF (BC_3L .GT. 0 .OR. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ) = 0.
        IF (BC_3U .GT. 0 .OR. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U)) = 0.
     END IF
     
  END IF
  
  !======================================================================================================================
  
  
  END SUBROUTINE exchange
  
  
  
  
  
  
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  SUBROUTINE exchange2(dir,vel_dir,SS1,SS2,SS3,NN1,NN2,NN3,phi) ! Anmerkung: Routine exakt identisch zu "exchange", allerdings mit variablen Intervallgrenzen
  
  IMPLICIT NONE
  
  INTEGER, INTENT(in   ) ::  dir
  INTEGER, INTENT(in   ) ::  vel_dir
  INTEGER, INTENT(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  ghost1LR( 1:b1U,SS2:NN2,SS3:NN3)
  REAL                   ::  ghost1LS( 1:b1U,SS2:NN2,SS3:NN3)
  REAL                   ::  ghost1UR(b1L:-1,SS2:NN2,SS3:NN3)
  REAL                   ::  ghost1US(b1L:-1,SS2:NN2,SS3:NN3)
  
  REAL                   ::  ghost2LR(SS1:NN1, 1:b2U,SS3:NN3)
  REAL                   ::  ghost2LS(SS1:NN1, 1:b2U,SS3:NN3)
  REAL                   ::  ghost2UR(SS1:NN1,b2L:-1,SS3:NN3)
  REAL                   ::  ghost2US(SS1:NN1,b2L:-1,SS3:NN3)
  
  REAL                   ::  ghost3LR(SS1:NN1,SS2:NN2, 1:b3U)
  REAL                   ::  ghost3LS(SS1:NN1,SS2:NN2, 1:b3U)
  REAL                   ::  ghost3UR(SS1:NN1,SS2:NN2,b3L:-1)
  REAL                   ::  ghost3US(SS1:NN1,SS2:NN2,b3L:-1)
  
  INTEGER                ::  length1L, length1U
  INTEGER                ::  length2L, length2U
  INTEGER                ::  length3L, length3U
  
  INTEGER                ::  i, j, k
  INTEGER                ::  dummy
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  IF (dir == 1) THEN
     
     length1L = (NN2-SS2+1)*(NN3-SS3+1)*ABS(b1L)
     length1U = (NN2-SS2+1)*(NN3-SS3+1)*ABS(b1U)
     
     IF (BC_1L == 0) CALL MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
     IF (BC_1U == 0) CALL MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
     
     IF (BC_1U == 0) THEN
        DO i = b1L, -1
           ghost1US(i,SS2:NN2,SS3:NN3) = phi((N1p+1+i),SS2:NN2,SS3:NN3)
        END DO
     END IF
     IF (BC_1L == 0) THEN
        DO i = 1, b1U
           ghost1LS(i,SS2:NN2,SS3:NN3) = phi((S1p-1+i),SS2:NN2,SS3:NN3)
        END DO
     END IF
     
     IF (BC_1U == 0) CALL MPI_SEND(ghost1US,length1L,MPI_REAL8,rank1U,1,COMM_CART,merror)
     IF (BC_1L == 0) CALL MPI_SEND(ghost1LS,length1U,MPI_REAL8,rank1L,2,COMM_CART,merror)
     
     IF (BC_1L == 0) CALL MPI_WAIT(req1L,status,merror)
     IF (BC_1U == 0) CALL MPI_WAIT(req1U,status,merror)
     
     IF (BC_1L == 0) CALL pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     IF (BC_1U == 0) CALL pseudocall(ghost1LR)
     
     IF (BC_1L == 0) THEN
        DO i = b1L, -1
           phi((S1p+i),SS2:NN2,SS3:NN3) = ghost1UR(i,SS2:NN2,SS3:NN3)
        END DO
     END IF
     IF (BC_1U == 0) THEN
        DO i = 1, b1U
           phi((N1p+i),SS2:NN2,SS3:NN3) = ghost1LR(i,SS2:NN2,SS3:NN3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (BC_1L == -1) THEN
        DO i = b1L, -1
           phi((S1p+i),SS2:NN2,SS3:NN3) = phi((N1p+1+i),SS2:NN2,SS3:NN3)
        END DO
        DO i = 1, b1U
           phi((N1p+i),SS2:NN2,SS3:NN3) = phi((S1p-1+i),SS2:NN2,SS3:NN3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (vel_dir == dir) THEN
        IF (BC_1L .GT. 0) phi( b1L  : -1     ,SS2:NN2,SS3:NN3) = 0.
        IF (BC_1U .GT. 0) phi((N1+1):(N1+b1U),SS2:NN2,SS3:NN3) = 0.
        
        IF (BC_1L  == -2) phi( b1L  :  0     ,SS2:NN2,SS3:NN3) = 0.
        IF (BC_1U  == -2) phi( N1   :(N1+b1U),SS2:NN2,SS3:NN3) = 0.
     ELSE
        IF (BC_1L .GT. 0 .OR. BC_1L == -2) phi( b1L  : 0      ,SS2:NN2,SS3:NN3) = 0.
        IF (BC_1U .GT. 0 .OR. BC_1U == -2) phi((N1+1):(N1+b1U),SS2:NN2,SS3:NN3) = 0.
     END IF
     
  END IF
  
  !======================================================================================================================
  
  IF (dir == 2) THEN
     
     length2L = (NN1-SS1+1)*(NN3-SS3+1)*ABS(b2L)
     length2U = (NN1-SS1+1)*(NN3-SS3+1)*ABS(b2U)
     
     IF (BC_2L == 0) CALL MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
     IF (BC_2U == 0) CALL MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
     
     IF (BC_2U == 0) THEN
        DO j = b2L, -1
           ghost2US(SS1:NN1,j,SS3:NN3) = phi(SS1:NN1,(N2p+1+j),SS3:NN3)
        END DO
     END IF
     IF (BC_2L == 0) THEN
        DO j = 1, b2U
           ghost2LS(SS1:NN1,j,SS3:NN3) = phi(SS1:NN1,(S2p-1+j),SS3:NN3)
        END DO
     END IF
     
     IF (BC_2U == 0) CALL MPI_SEND(ghost2US,length2L,MPI_REAL8,rank2U,3,COMM_CART,merror)
     IF (BC_2L == 0) CALL MPI_SEND(ghost2LS,length2U,MPI_REAL8,rank2L,4,COMM_CART,merror)
     
     IF (BC_2L == 0) CALL MPI_WAIT(req2L,status,merror)
     IF (BC_2U == 0) CALL MPI_WAIT(req2U,status,merror)
     
     IF (BC_2L == 0) CALL pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     IF (BC_2U == 0) CALL pseudocall(ghost2LR)
     
     IF (BC_2L == 0) THEN
        DO j = b2L, -1
           phi(SS1:NN1,(S2p+j),SS3:NN3) = ghost2UR(SS1:NN1,j,SS3:NN3)
        END DO
     END IF
     IF (BC_2U == 0) THEN
        DO j = 1, b2U
           phi(SS1:NN1,(N2p+j),SS3:NN3) = ghost2LR(SS1:NN1,j,SS3:NN3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (BC_2L == -1) THEN
        DO j = b2L, -1
           phi(SS1:NN1,(S2p+j),SS3:NN3) = phi(SS1:NN1,(N2p+1+j),SS3:NN3)
        END DO
        DO j = 1, b2U
           phi(SS1:NN1,(N2p+j),SS3:NN3) = phi(SS1:NN1,(S2p-1+j),SS3:NN3)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (vel_dir == dir) THEN
        IF (BC_2L .GT. 0) phi(SS1:NN1, b2L  : -1     ,SS3:NN3) = 0.
        IF (BC_2U .GT. 0) phi(SS1:NN1,(N2+1):(N2+b2U),SS3:NN3) = 0.
        
        IF (BC_2L  == -2) phi(SS1:NN1, b2L  :  0     ,SS3:NN3) = 0.
        IF (BC_2U  == -2) phi(SS1:NN1, N2   :(N2+b2U),SS3:NN3) = 0.
     ELSE
        IF (BC_2L .GT. 0 .OR. BC_2L == -2) phi(SS1:NN1, b2L  : 0      ,SS3:NN3) = 0.
        IF (BC_2U .GT. 0 .OR. BC_2U == -2) phi(SS1:NN1,(N2+1):(N2+b2U),SS3:NN3) = 0.
     END IF
     
  END IF
  
  !======================================================================================================================
  
  IF (dir == 3 .AND. dimens == 3) THEN
     
     length3L = (NN1-SS1+1)*(NN2-SS2+1)*ABS(b3L)
     length3U = (NN1-SS1+1)*(NN2-SS2+1)*ABS(b3U)
     
     IF (BC_3L == 0) CALL MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
     IF (BC_3U == 0) CALL MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
     
     IF (BC_3U == 0) THEN
        DO k = b3L, -1
           ghost3US(SS1:NN1,SS2:NN2,k) = phi(SS1:NN1,SS2:NN2,(N3p+1+k))
        END DO
     END IF
     IF (BC_3L == 0) THEN
        DO k = 1, b3U
           ghost3LS(SS1:NN1,SS2:NN2,k) = phi(SS1:NN1,SS2:NN2,(S3p-1+k))
        END DO
     END IF
     
     IF (BC_3U == 0) CALL MPI_SEND(ghost3US,length3L,MPI_REAL8,rank3U,5,COMM_CART,merror)
     IF (BC_3L == 0) CALL MPI_SEND(ghost3LS,length3U,MPI_REAL8,rank3L,6,COMM_CART,merror)
     
     IF (BC_3L == 0) CALL MPI_WAIT(req3L,status,merror)
     IF (BC_3U == 0) CALL MPI_WAIT(req3U,status,merror)
     
     IF (BC_3L == 0) CALL pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
     IF (BC_3U == 0) CALL pseudocall(ghost3LR)
     
     IF (BC_3L == 0) THEN
        DO k = b3L, -1
           phi(SS1:NN1,SS2:NN2,(S3p+k)) = ghost3UR(SS1:NN1,SS2:NN2,k)
        END DO
     END IF
     IF (BC_3U == 0) THEN
        DO k = 1, b3U
           phi(SS1:NN1,SS2:NN2,(N3p+k)) = ghost3LR(SS1:NN1,SS2:NN2,k)
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (BC_3L == -1) THEN
        DO k = b3L, -1
           phi(SS1:NN1,SS2:NN2,(S3p+k)) = phi(SS1:NN1,SS2:NN2,(N3p+1+k))
        END DO
        DO k = 1, b3U
           phi(SS1:NN1,SS2:NN2,(N3p+k)) = phi(SS1:NN1,SS2:NN2,(S3p-1+k))
        END DO
     END IF
     
     !-------------------------------------------------------------------------------------------------------------------
     
     IF (vel_dir == dir) THEN
        IF (BC_3L .GT. 0) phi(SS1:NN1,SS2:NN2, b3L  : -1     ) = 0.
        IF (BC_3U .GT. 0) phi(SS1:NN1,SS2:NN2,(N3+1):(N3+b3U)) = 0.
        
        IF (BC_3L  == -2) phi(SS1:NN1,SS2:NN2, b3L  :  0     ) = 0.
        IF (BC_3U  == -2) phi(SS1:NN1,SS2:NN2, N3   :(N3+b3U)) = 0.
     ELSE
        IF (BC_3L .GT. 0 .OR. BC_3L == -2) phi(SS1:NN1,SS2:NN2, b3L  : 0      ) = 0.
        IF (BC_3U .GT. 0 .OR. BC_3U == -2) phi(SS1:NN1,SS2:NN2,(N3+1):(N3+b3U)) = 0.
     END IF
     
  END IF
  
  !======================================================================================================================
  
  
  END SUBROUTINE exchange2
  
  
  
  
  
  
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  SUBROUTINE exchange_all_all(vel_yes,phi)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(in)    ::  vel_yes
  
  REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  REAL                   ::  ghost1LR( 1:b1U,1:N2,1:N3,1:3)
  REAL                   ::  ghost1LS( 1:b1U,1:N2,1:N3,1:3)
  REAL                   ::  ghost1UR(b1L:-1,1:N2,1:N3,1:3)
  REAL                   ::  ghost1US(b1L:-1,1:N2,1:N3,1:3)
  
  REAL                   ::  ghost2LR(1:N1, 1:b2U,1:N3,1:3)
  REAL                   ::  ghost2LS(1:N1, 1:b2U,1:N3,1:3)
  REAL                   ::  ghost2UR(1:N1,b2L:-1,1:N3,1:3)
  REAL                   ::  ghost2US(1:N1,b2L:-1,1:N3,1:3)
  
  REAL                   ::  ghost3LR(1:N1,1:N2, 1:b3U,1:3)
  REAL                   ::  ghost3LS(1:N1,1:N2, 1:b3U,1:3)
  REAL                   ::  ghost3UR(1:N1,1:N2,b3L:-1,1:3)
  REAL                   ::  ghost3US(1:N1,1:N2,b3L:-1,1:3)
  
  INTEGER                ::  length1L, length1U
  INTEGER                ::  length2L, length2U
  INTEGER                ::  length3L, length3U
  
  INTEGER                ::  i, j, k
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
  !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
  !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
  !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
  !                kein Prefetch einbaut.                                                                               !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  length1L = N2*N3*ABS(b1L)*3
  length1U = N2*N3*ABS(b1U)*3
  
  IF (BC_1L == 0) CALL MPI_IRECV(ghost1UR,length1L,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
  IF (BC_1U == 0) CALL MPI_IRECV(ghost1LR,length1U,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
  
  IF (BC_1U == 0) THEN
     DO i = b1L, -1
        ghost1US(i,1:N2,1:N3,:) = phi((N1p+1+i),1:N2,1:N3,:) ! TEST!!! Hier koennte man auch 1:dimens schreiben!!!
     END DO
  END IF
  IF (BC_1L == 0) THEN
     DO i = 1, b1U
        ghost1LS(i,1:N2,1:N3,:) = phi((S1p-1+i),1:N2,1:N3,:)
     END DO
  END IF
  
  IF (BC_1U == 0) CALL MPI_SEND(ghost1US,length1L,MPI_REAL8,rank1U,1,COMM_CART,merror)
  IF (BC_1L == 0) CALL MPI_SEND(ghost1LS,length1U,MPI_REAL8,rank1L,2,COMM_CART,merror)
  
  IF (BC_1L == 0) CALL MPI_WAIT(req1L,status,merror)
  IF (BC_1U == 0) CALL MPI_WAIT(req1U,status,merror)
  
  IF (BC_1L == 0) CALL pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (BC_1U == 0) CALL pseudocall(ghost1LR)
  
  IF (BC_1L == 0) THEN
     DO i = b1L, -1
        phi((S1p+i),1:N2,1:N3,:) = ghost1UR(i,1:N2,1:N3,:)
     END DO
  END IF
  IF (BC_1U == 0) THEN
     DO i = 1, b1U
        phi((N1p+i),1:N2,1:N3,:) = ghost1LR(i,1:N2,1:N3,:)
     END DO
  END IF
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (BC_1L == -1) THEN
     DO i = b1L, -1
        phi((S1p+i),1:N2,1:N3,:) = phi((N1p+1+i),1:N2,1:N3,:)
     END DO
     DO i = 1, b1U
        phi((N1p+i),1:N2,1:N3,:) = phi((S1p-1+i),1:N2,1:N3,:)
     END DO
  END IF
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (vel_yes) THEN
     IF (BC_1L .GT. 0) phi( b1L  : -1     ,1:N2,1:N3,1) = 0.
     IF (BC_1U .GT. 0) phi((N1+1):(N1+b1U),1:N2,1:N3,1) = 0.
     
     IF (BC_1L  == -2) phi( b1L  :  0     ,1:N2,1:N3,1) = 0.
     IF (BC_1U  == -2) phi( N1   :(N1+b1U),1:N2,1:N3,1) = 0.
  ELSE
     IF (BC_1L .GT. 0 .OR. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3,1) = 0.
     IF (BC_1U .GT. 0 .OR. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3,1) = 0.
  END IF
  
  IF (BC_1L .GT. 0 .OR. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3,2) = 0.
  IF (BC_1L .GT. 0 .OR. BC_1L == -2) phi( b1L  : 0      ,1:N2,1:N3,3) = 0.
  IF (BC_1U .GT. 0 .OR. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3,2) = 0.
  IF (BC_1U .GT. 0 .OR. BC_1U == -2) phi((N1+1):(N1+b1U),1:N2,1:N3,3) = 0.
  
  !======================================================================================================================
  
  length2L = N1*N3*ABS(b2L)*3
  length2U = N1*N3*ABS(b2U)*3
  
  IF (BC_2L == 0) CALL MPI_IRECV(ghost2UR,length2L,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
  IF (BC_2U == 0) CALL MPI_IRECV(ghost2LR,length2U,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
  
  IF (BC_2U == 0) THEN
     DO j = b2L, -1
        ghost2US(1:N1,j,1:N3,:) = phi(1:N1,(N2p+1+j),1:N3,:)
     END DO
  END IF
  IF (BC_2L == 0) THEN
     DO j = 1, b2U
        ghost2LS(1:N1,j,1:N3,:) = phi(1:N1,(S2p-1+j),1:N3,:)
     END DO
  END IF
  
  IF (BC_2U == 0) CALL MPI_SEND(ghost2US,length2L,MPI_REAL8,rank2U,3,COMM_CART,merror)
  IF (BC_2L == 0) CALL MPI_SEND(ghost2LS,length2U,MPI_REAL8,rank2L,4,COMM_CART,merror)
  
  IF (BC_2L == 0) CALL MPI_WAIT(req2L,status,merror)
  IF (BC_2U == 0) CALL MPI_WAIT(req2U,status,merror)
  
  IF (BC_2L == 0) CALL pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (BC_2U == 0) CALL pseudocall(ghost2LR)
  
  IF (BC_2L == 0) THEN
     DO j = b2L, -1
        phi(1:N1,(S2p+j),1:N3,:) = ghost2UR(1:N1,j,1:N3,:)
     END DO
  END IF
  IF (BC_2U == 0) THEN
     DO j = 1, b2U
        phi(1:N1,(N2p+j),1:N3,:) = ghost2LR(1:N1,j,1:N3,:)
     END DO
  END IF
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (BC_2L == -1) THEN
     DO j = b2L, -1
        phi(1:N1,(S2p+j),1:N3,:) = phi(1:N1,(N2p+1+j),1:N3,:)
     END DO
     DO j = 1, b2U
        phi(1:N1,(N2p+j),1:N3,:) = phi(1:N1,(S2p-1+j),1:N3,:)
     END DO
  END IF
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (vel_yes) THEN
     IF (BC_2L .GT. 0) phi(1:N1, b2L  : -1     ,1:N3,2) = 0.
     IF (BC_2U .GT. 0) phi(1:N1,(N2+1):(N2+b2U),1:N3,2) = 0.
     
     IF (BC_2L  == -2) phi(1:N1, b2L  :  0     ,1:N3,2) = 0.
     IF (BC_2U  == -2) phi(1:N1, N2   :(N2+b2U),1:N3,2) = 0.
  ELSE
     IF (BC_2L .GT. 0 .OR. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3,2) = 0.
     IF (BC_2U .GT. 0 .OR. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3,2) = 0.
  END IF
  
  IF (BC_2L .GT. 0 .OR. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3,1) = 0.
  IF (BC_2L .GT. 0 .OR. BC_2L == -2) phi(1:N1, b2L  : 0      ,1:N3,3) = 0.
  IF (BC_2U .GT. 0 .OR. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3,1) = 0.
  IF (BC_2U .GT. 0 .OR. BC_2U == -2) phi(1:N1,(N2+1):(N2+b2U),1:N3,3) = 0.
  
  !======================================================================================================================
  
  IF (dimens == 3) THEN
  
  length3L = N1*N2*ABS(b3L)*3
  length3U = N1*N2*ABS(b3U)*3
  
  IF (BC_3L == 0) CALL MPI_IRECV(ghost3UR,length3L,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
  IF (BC_3U == 0) CALL MPI_IRECV(ghost3LR,length3U,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
  
  IF (BC_3U == 0) THEN
     DO k = b3L, -1
        ghost3US(1:N1,1:N2,k,:) = phi(1:N1,1:N2,(N3p+1+k),:)
     END DO
  END IF
  IF (BC_3L == 0) THEN
     DO k = 1, b3U
        ghost3LS(1:N1,1:N2,k,:) = phi(1:N1,1:N2,(S3p-1+k),:)
     END DO
  END IF
  
  IF (BC_3U == 0) CALL MPI_SEND(ghost3US,length3L,MPI_REAL8,rank3U,5,COMM_CART,merror)
  IF (BC_3L == 0) CALL MPI_SEND(ghost3LS,length3U,MPI_REAL8,rank3L,6,COMM_CART,merror)
  
  IF (BC_3L == 0) CALL MPI_WAIT(req3L,status,merror)
  IF (BC_3U == 0) CALL MPI_WAIT(req3U,status,merror)
  
  IF (BC_3L == 0) CALL pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (BC_3U == 0) CALL pseudocall(ghost3LR)
  
  IF (BC_3L == 0) THEN
     DO k = b3L, -1
        phi(1:N1,1:N2,(S3p+k),:) = ghost3UR(1:N1,1:N2,k,:)
     END DO
  END IF
  IF (BC_3U == 0) THEN
     DO k = 1, b3U
        phi(1:N1,1:N2,(N3p+k),:) = ghost3LR(1:N1,1:N2,k,:)
     END DO
  END IF
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (BC_3L == -1) THEN
     DO k = b3L, -1
        phi(1:N1,1:N2,(S3p+k),:) = phi(1:N1,1:N2,(N3p+1+k),:)
     END DO
     DO k = 1, b3U
        phi(1:N1,1:N2,(N3p+k),:) = phi(1:N1,1:N2,(S3p-1+k),:)
     END DO
  END IF
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (vel_yes) THEN
     IF (BC_3L .GT. 0) phi(1:N1,1:N2, b3L  : -1     ,3) = 0.
     IF (BC_3U .GT. 0) phi(1:N1,1:N2,(N3+1):(N3+b3U),3) = 0.
     
     IF (BC_3L  == -2) phi(1:N1,1:N2, b3L  :  0     ,3) = 0.
     IF (BC_3U  == -2) phi(1:N1,1:N2, N3   :(N3+b3U),3) = 0.
  ELSE
     IF (BC_3L .GT. 0 .OR. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ,3) = 0.
     IF (BC_3U .GT. 0 .OR. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U),3) = 0.
  END IF
  
  IF (BC_3L .GT. 0 .OR. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ,1) = 0.
  IF (BC_3L .GT. 0 .OR. BC_3L == -2) phi(1:N1,1:N2, b3L  : 0      ,2) = 0.
  IF (BC_3U .GT. 0 .OR. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U),1) = 0.
  IF (BC_3U .GT. 0 .OR. BC_3U == -2) phi(1:N1,1:N2,(N3+1):(N3+b3U),2) = 0.
  
  END IF
  !======================================================================================================================
  
  
  END SUBROUTINE exchange_all_all
  
  
  
  
  
  
  
  
  
!pgi$r nodepchk
!pgi$r ivdep
  SUBROUTINE exchange_relax(g,fb1,fb2,fb3,vel_dir,mirror_yes,phi) ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  !*************************************************************************************************
  INTEGER, INTENT(in   ) ::  g
  REAL   , INTENT(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  INTEGER                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  INTEGER                ::  rank1L, rank1U, rank2L, rank2U, rank3L, rank3U
  !*************************************************************************************************
  
  INTEGER, INTENT(in)    ::  fb1
  INTEGER, INTENT(in)    ::  fb2
  INTEGER, INTENT(in)    ::  fb3
  
  INTEGER, INTENT(in)    ::  vel_dir
  LOGICAL, INTENT(in)    ::  mirror_yes
  
  REAL                   ::  ghost1LR(1:NN(2,g),1:NN(3,g)) ! TEST!!! Arrays sollten frueher oder spaeter global allociert werden ...
  REAL                   ::  ghost1LS(1:NN(2,g),1:NN(3,g))
  REAL                   ::  ghost1UR(1:NN(2,g),1:NN(3,g))
  REAL                   ::  ghost1US(1:NN(2,g),1:NN(3,g))
  
  REAL                   ::  ghost2LR(1:NN(1,g),1:NN(3,g))
  REAL                   ::  ghost2LS(1:NN(1,g),1:NN(3,g))
  REAL                   ::  ghost2UR(1:NN(1,g),1:NN(3,g))
  REAL                   ::  ghost2US(1:NN(1,g),1:NN(3,g))
  
  REAL                   ::  ghost3LR(1:NN(1,g),1:NN(2,g))
  REAL                   ::  ghost3LS(1:NN(1,g),1:NN(2,g))
  REAL                   ::  ghost3UR(1:NN(1,g),1:NN(2,g))
  REAL                   ::  ghost3US(1:NN(1,g),1:NN(2,g))
  
  INTEGER                ::  length1
  INTEGER                ::  length2
  INTEGER                ::  length3
  
  INTEGER                ::  i, N1, S1R, N1R
  INTEGER                ::  j, N2, S2R, N2R
  INTEGER                ::  k, N3, S3R, N3R
  
  
  !*****************************************************************************************
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  BC_1L = BC(1,1,g)
  BC_1U = BC(2,1,g)
  BC_2L = BC(1,2,g)
  BC_2U = BC(2,2,g)
  BC_3L = BC(1,3,g)
  BC_3U = BC(2,3,g)
  
  rank1L = ngb(1,1,g)
  rank1U = ngb(2,1,g)
  rank2L = ngb(1,2,g)
  rank2U = ngb(2,2,g)
  rank3L = ngb(1,3,g)
  rank3U = ngb(2,3,g)
  
  S1R  = SNB(1,1,g)
  S2R  = SNB(1,2,g)
  S3R  = SNB(1,3,g)
  
  N1R  = SNB(2,1,g)
  N2R  = SNB(2,2,g)
  N3R  = SNB(2,3,g)
  !*****************************************************************************************
  
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Intervallgrenzen sind der Einfachheit halber analog zu exchange gewählt.                             !
  !              - Im Fall von Symmetrie-RB muss IMMER entsprechend oben und/oder unten gespiegelt werden.              !
  !              - Die Spiegelung wird anstelle einer Koeffizienten-Manipulation vorgenommen, da ansonsten auch in den  !
  !                Grob-Gitter-Routinen die Punkte ausserhalb der betrachteten Domain zu Null gesetzt werden müssten,   !
  !                was im Endeffekt (vermutlich) der gleiche Aufwand wäre.                                              !
  !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
  !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
  !              - Die "vel_dir"-Variante wird ausschliesslich von allen Helmholtz-Operationen auf dem feinsten Gitter  !
  !                benutzt, INKL. Interpolation.                                                                        !
  !              - Bei Spiegelung (mirror_yes bzw. vel_dir) wird das Nullsetzen des Randes einer expliziten Behandlung  !
  !                in den aufrufenden Routinen vorgezogen, da die Komplexität und der Programmieraufwand im Vergleich   !
  !                zum Geschwindigkeitsgewinn zu gross erscheinen (im Gegensatz zum Poisson-Problem mit Neumann-RB, bei !
  !                dem am Rand kein Laplace-3d-Stencil sondern nur eine 1d-Ableitung vorkommt (!)).                     !
  !              - Die "mirror_yes"-Variante wird ausschliesslich von allen Helmholtz- und Poisson-Operationen benutzt, !
  !                OHNE Interpolation.                                                                                  !
  !              - In der ersten Richtungen kann der Compiler offenbar nicht vektorisieren und in der zweiten zudem     !
  !                nicht unrollen. Der Grund dafür ist nicht klar.                                                      !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  !======================================================================================================================
  
  length1 = N2*N3
  
  IF (fb1 .GE. 0 .AND. BC_1L ==  0) CALL MPI_IRECV(ghost1UR,length1,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
  IF (fb1 .LE. 0 .AND. BC_1U ==  0) CALL MPI_IRECV(ghost1LR,length1,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
  
  IF (fb1 .GE. 0 .AND. BC_1U ==  0) ghost1US(1:N2,1:N3) = phi(N1R,1:N2,1:N3)
  IF (fb1 .LE. 0 .AND. BC_1L ==  0) ghost1LS(1:N2,1:N3) = phi(S1R,1:N2,1:N3)
  
  IF (fb1 .GE. 0 .AND. BC_1U ==  0) CALL MPI_SEND(ghost1US,length1,MPI_REAL8,rank1U,1,COMM_CART,merror)
  IF (fb1 .LE. 0 .AND. BC_1L ==  0) CALL MPI_SEND(ghost1LS,length1,MPI_REAL8,rank1L,2,COMM_CART,merror)
  
  IF (fb1 .GE. 0 .AND. BC_1L ==  0) CALL MPI_WAIT(req1L,status,merror)
  IF (fb1 .LE. 0 .AND. BC_1U ==  0) CALL MPI_WAIT(req1U,status,merror)
  
  IF (fb1 .GE. 0 .AND. BC_1L ==  0) CALL pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (fb1 .LE. 0 .AND. BC_1U ==  0) CALL pseudocall(ghost1LR)
  
  IF (fb1 .GE. 0 .AND. BC_1L ==  0) phi(S1R-1,1:N2,1:N3) = ghost1UR(1:N2,1:N3)
  IF (fb1 .LE. 0 .AND. BC_1U ==  0) phi(N1R+1,1:N2,1:N3) = ghost1LR(1:N2,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (fb1 .GE. 0 .AND. BC_1L == -1) phi(S1R-1,1:N2,1:N3) = phi(N1R,1:N2,1:N3)
  IF (fb1 .LE. 0 .AND. BC_1L == -1) phi(N1R+1,1:N2,1:N3) = phi(S1R,1:N2,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (vel_dir == 1) THEN
     IF (BC_1L == -2) phi(0   ,1:N2,1:N3) = 0.
     IF (BC_1U == -2) phi(N1  ,1:N2,1:N3) = 0.
  ELSE IF (mirror_yes) THEN
     IF (BC_1L == -2) phi(0   ,1:N2,1:N3) = 0.
     IF (BC_1U == -2) phi(N1+1,1:N2,1:N3) = 0.
  END IF
  
  !======================================================================================================================
  
  length2 = N1*N3
  
  IF (fb2 .GE. 0 .AND. BC_2L ==  0) CALL MPI_IRECV(ghost2UR,length2,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
  IF (fb2 .LE. 0 .AND. BC_2U ==  0) CALL MPI_IRECV(ghost2LR,length2,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
  
  IF (fb2 .GE. 0 .AND. BC_2U ==  0) ghost2US(1:N1,1:N3) = phi(1:N1,N2R,1:N3)
  IF (fb2 .LE. 0 .AND. BC_2L ==  0) ghost2LS(1:N1,1:N3) = phi(1:N1,S2R,1:N3)
  
  IF (fb2 .GE. 0 .AND. BC_2U ==  0) CALL MPI_SEND(ghost2US,length2,MPI_REAL8,rank2U,3,COMM_CART,merror)
  IF (fb2 .LE. 0 .AND. BC_2L ==  0) CALL MPI_SEND(ghost2LS,length2,MPI_REAL8,rank2L,4,COMM_CART,merror)
  
  IF (fb2 .GE. 0 .AND. BC_2L ==  0) CALL MPI_WAIT(req2L,status,merror)
  IF (fb2 .LE. 0 .AND. BC_2U ==  0) CALL MPI_WAIT(req2U,status,merror)
  
  IF (fb2 .GE. 0 .AND. BC_2L ==  0) CALL pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (fb2 .LE. 0 .AND. BC_2U ==  0) CALL pseudocall(ghost2LR)
  
  IF (fb2 .GE. 0 .AND. BC_2L ==  0) phi(1:N1,S2R-1,1:N3) = ghost2UR(1:N1,1:N3)
  IF (fb2 .LE. 0 .AND. BC_2U ==  0) phi(1:N1,N2R+1,1:N3) = ghost2LR(1:N1,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (fb2 .GE. 0 .AND. BC_2L == -1) phi(1:N1,S2R-1,1:N3) = phi(1:N1,N2R,1:N3)
  IF (fb2 .LE. 0 .AND. BC_2L == -1) phi(1:N1,N2R+1,1:N3) = phi(1:N1,S2R,1:N3)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (vel_dir == 2) THEN
     IF (BC_2L == -2) phi(1:N1,0   ,1:N3) = 0.
     IF (BC_2U == -2) phi(1:N1,N2  ,1:N3) = 0.
  ELSE IF (mirror_yes) THEN
     IF (BC_2L == -2) phi(1:N1,0   ,1:N3) = 0.
     IF (BC_2U == -2) phi(1:N1,N2+1,1:N3) = 0.
  END IF
  
  !======================================================================================================================
  IF (dimens == 3) THEN
  
  length3 = N1*N2
  
  IF (fb3 .GE. 0 .AND. BC_3L ==  0) CALL MPI_IRECV(ghost3UR,length3,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
  IF (fb3 .LE. 0 .AND. BC_3U ==  0) CALL MPI_IRECV(ghost3LR,length3,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
  
  IF (fb3 .GE. 0 .AND. BC_3U ==  0) ghost3US(1:N1,1:N2) = phi(1:N1,1:N2,N3R)
  IF (fb3 .LE. 0 .AND. BC_3L ==  0) ghost3LS(1:N1,1:N2) = phi(1:N1,1:N2,S3R)
  
  IF (fb3 .GE. 0 .AND. BC_3U ==  0) CALL MPI_SEND(ghost3US,length3,MPI_REAL8,rank3U,5,COMM_CART,merror)
  IF (fb3 .LE. 0 .AND. BC_3L ==  0) CALL MPI_SEND(ghost3LS,length3,MPI_REAL8,rank3L,6,COMM_CART,merror)
  
  IF (fb3 .GE. 0 .AND. BC_3L ==  0) CALL MPI_WAIT(req3L,status,merror)
  IF (fb3 .LE. 0 .AND. BC_3U ==  0) CALL MPI_WAIT(req3U,status,merror)
  
  IF (fb3 .GE. 0 .AND. BC_3L ==  0) CALL pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
  IF (fb3 .LE. 0 .AND. BC_3U ==  0) CALL pseudocall(ghost3LR)
  
  IF (fb3 .GE. 0 .AND. BC_3L ==  0) phi(1:N1,1:N2,S3R-1) = ghost3UR(1:N1,1:N2)
  IF (fb3 .LE. 0 .AND. BC_3U ==  0) phi(1:N1,1:N2,N3R+1) = ghost3LR(1:N1,1:N2)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (fb3 .GE. 0 .AND. BC_3L == -1) phi(1:N1,1:N2,S3R-1) = phi(1:N1,1:N2,N3R)
  IF (fb3 .LE. 0 .AND. BC_3L == -1) phi(1:N1,1:N2,N3R+1) = phi(1:N1,1:N2,S3R)
  
  !----------------------------------------------------------------------------------------------------------------------
  
  IF (vel_dir == 3) THEN
     IF (BC_3L == -2) phi(1:N1,1:N2,0   ) = 0.
     IF (BC_3U == -2) phi(1:N1,1:N2,N3  ) = 0.
  ELSE IF (mirror_yes) THEN
     IF (BC_3L == -2) phi(1:N1,1:N2,0   ) = 0.
     IF (BC_3U == -2) phi(1:N1,1:N2,N3+1) = 0.
  END IF
  
  END IF
  !======================================================================================================================
  
  
  END SUBROUTINE exchange_relax
  
  
  
  
  
  
  
  
  
  
  
   SUBROUTINE exchange_part(dir,vel_dir,ss,phi)
   
   IMPLICIT NONE
   
   INTEGER, INTENT(in   ) ::  dir
   INTEGER, INTENT(in   ) ::  vel_dir
   INTEGER, INTENT(in   ) ::  ss
   
   REAL   , INTENT(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
   
   REAL                   ::  ghost1LR(0:(N2+1),0:(N3+1))
   REAL                   ::  ghost1LS(0:(N2+1),0:(N3+1))
   REAL                   ::  ghost1UR(0:(N2+1),0:(N3+1))
   REAL                   ::  ghost1US(0:(N2+1),0:(N3+1))
   
   REAL                   ::  ghost2LR(0:(N1+1),0:(N3+1))
   REAL                   ::  ghost2LS(0:(N1+1),0:(N3+1))
   REAL                   ::  ghost2UR(0:(N1+1),0:(N3+1))
   REAL                   ::  ghost2US(0:(N1+1),0:(N3+1))
   
   REAL                   ::  ghost3LR(0:(N1+1),0:(N2+1))
   REAL                   ::  ghost3LS(0:(N1+1),0:(N2+1))
   REAL                   ::  ghost3UR(0:(N1+1),0:(N2+1))
   REAL                   ::  ghost3US(0:(N1+1),0:(N2+1))
   
   INTEGER                ::  length
   INTEGER                ::  i, j, k
   
   
   !---------------------------------------------------------------------------------------------------------------------!
   ! Anmerkungen: - [Sij,Nij] ist Untermenge von [Sip,Nip]                                                               !
   !              - [Sip,Nip] ist Untermenge von [1,Ni]                                                                  !
   !              - [1,Ni] hat den Vorteil, dass Intervallgrenzen fest einprogrammiert sind,                             !
   !                was prinzipiell zu einem Vorteil bei der Effizienz führen sollte.                                    !
   !              - Bei Spiegelung muss die zur Spiegelungsebene orthogonale Geschwindigkeitskomponente auch im          !
   !                Vorzeichen gespiegelt werden, da die Spiegelungsebene nicht durchströmt werden kann!.                !
   !              - Alle Kopier-Schleifen sind ausgeschrieben, da der PGI-Compiler ansonsten keine Vektorisierung bzw.   !
   !                kein Prefetch einbaut.                                                                               !
   !---------------------------------------------------------------------------------------------------------------------!
   
   
   !======================================================================================================================
   
   IF (dir == 1) THEN
      
      length = (N2+2)*(N3+2)
      
      IF (BC_1L == 0) CALL MPI_IRECV(ghost1UR,length,MPI_REAL8,rank1L,1,COMM_CART,req1L,merror)
      IF (BC_1U == 0) CALL MPI_IRECV(ghost1LR,length,MPI_REAL8,rank1U,2,COMM_CART,req1U,merror)
      
      IF (BC_1U == 0) ghost1US(0:(N2+1),0:(N3+1)) = phi(N1p+ss,0:(N2+1),0:(N3+1))
      IF (BC_1L == 0) ghost1LS(0:(N2+1),0:(N3+1)) = phi(S1p-ss,0:(N2+1),0:(N3+1))
      
      IF (BC_1U == 0) CALL MPI_SEND(ghost1US,length,MPI_REAL8,rank1U,1,COMM_CART,merror)
      IF (BC_1L == 0) CALL MPI_SEND(ghost1LS,length,MPI_REAL8,rank1L,2,COMM_CART,merror)
      
      IF (BC_1L == 0) CALL MPI_WAIT(req1L,status,merror)
      IF (BC_1U == 0) CALL MPI_WAIT(req1U,status,merror)
      
      IF (BC_1L == 0) CALL pseudocall(ghost1UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
      IF (BC_1U == 0) CALL pseudocall(ghost1LR)
      
      IF (ss == 1) THEN
         IF (BC_1L == 0) phi(S1p  ,0:(N2+1),0:(N3+1)) = ghost1UR(0:(N2+1),0:(N3+1)) + phi(S1p,0:(N2+1),0:(N3+1))
         IF (BC_1U == 0) phi(N1p  ,0:(N2+1),0:(N3+1)) = ghost1LR(0:(N2+1),0:(N3+1)) + phi(N1p,0:(N2+1),0:(N3+1))
      ELSE
         IF (BC_1L == 0) phi(S1p-1,0:(N2+1),0:(N3+1)) = ghost1UR(0:(N2+1),0:(N3+1))
         IF (BC_1U == 0) phi(N1p+1,0:(N2+1),0:(N3+1)) = ghost1LR(0:(N2+1),0:(N3+1))
      END IF
      
      !-------------------------------------------------------------------------------------------------------------------
      
      IF (ss == 1) THEN
         IF (BC_1L == -1) THEN
            phi(S1p  ,0:(N2+1),0:(N3+1)) = phi(S1p,0:(N2+1),0:(N3+1)) + phi(N1p+1,0:(N2+1),0:(N3+1))
            phi(N1p  ,0:(N2+1),0:(N3+1)) = phi(N1p,0:(N2+1),0:(N3+1)) + phi(S1p-1,0:(N2+1),0:(N3+1))
         END IF
      ELSE
         IF (BC_1L == -1) THEN
            phi(S1p-1,0:(N2+1),0:(N3+1)) = phi(N1p,0:(N2+1),0:(N3+1))
            phi(N1p+1,0:(N2+1),0:(N3+1)) = phi(S1p,0:(N2+1),0:(N3+1))
         END IF
      END IF
      
      !-------------------------------------------------------------------------------------------------------------------
      
      IF (vel_dir == dir) THEN
         IF (BC_1L .GT. 0) phi(  -1,0:(N2+1),0:(N3+1)) = 0.
         IF (BC_1U .GT. 0) phi(N1+1,0:(N2+1),0:(N3+1)) = 0.
         
         IF (BC_1L  == -2) phi(   0,0:(N2+1),0:(N3+1)) = -phi(   1,0:(N2+1),0:(N3+1))
         IF (BC_1U  == -2) phi(N1  ,0:(N2+1),0:(N3+1)) = -phi(N1-1,0:(N2+1),0:(N3+1))
      ELSE
         IF (BC_1L .GT. 0) phi(   0,0:(N2+1),0:(N3+1)) = 0.
         IF (BC_1U .GT. 0) phi(N1+1,0:(N2+1),0:(N3+1)) = 0.
         
         !IF (ss == 1) THEN
         !   IF (BC_1L == -2) phi(1 ,0:(N2+1),0:(N3+1)) = 2.*phi(1 ,0:(N2+1),0:(N3+1)) ! TEST!!! feedback forces
         !   IF (BC_1U == -2) phi(N1,0:(N2+1),0:(N3+1)) = 2.*phi(N1,0:(N2+1),0:(N3+1))
         !ELSE
            IF (BC_1L == -2) phi(   0,0:(N2+1),0:(N3+1)) = 0. ! TEST!!! Fluid
            IF (BC_1U == -2) phi(N1+1,0:(N2+1),0:(N3+1)) = 0.
         !END IF
      END IF
      
   END IF
   
   !======================================================================================================================
   
   IF (dir == 2) THEN
      
      length = (N1+2)*(N3+2)
      
      IF (BC_2L == 0) CALL MPI_IRECV(ghost2UR,length,MPI_REAL8,rank2L,3,COMM_CART,req2L,merror)
      IF (BC_2U == 0) CALL MPI_IRECV(ghost2LR,length,MPI_REAL8,rank2U,4,COMM_CART,req2U,merror)
      
      IF (BC_2U == 0) ghost2US(0:(N1+1),0:(N3+1)) = phi(0:(N1+1),N2p+ss,0:(N3+1))
      IF (BC_2L == 0) ghost2LS(0:(N1+1),0:(N3+1)) = phi(0:(N1+1),S2p-ss,0:(N3+1))
      
      IF (BC_2U == 0) CALL MPI_SEND(ghost2US,length,MPI_REAL8,rank2U,3,COMM_CART,merror)
      IF (BC_2L == 0) CALL MPI_SEND(ghost2LS,length,MPI_REAL8,rank2L,4,COMM_CART,merror)
      
      IF (BC_2L == 0) CALL MPI_WAIT(req2L,status,merror)
      IF (BC_2U == 0) CALL MPI_WAIT(req2U,status,merror)
      
      IF (BC_2L == 0) CALL pseudocall(ghost2UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
      IF (BC_2U == 0) CALL pseudocall(ghost2LR)
      
      IF (ss == 1) THEN
         IF (BC_2L == 0) phi(0:(N1+1),S2p  ,0:(N3+1)) = ghost2UR(0:(N1+1),0:(N3+1)) + phi(0:(N1+1),S2p,0:(N3+1))
         IF (BC_2U == 0) phi(0:(N1+1),N2p  ,0:(N3+1)) = ghost2LR(0:(N1+1),0:(N3+1)) + phi(0:(N1+1),N2p,0:(N3+1))
      ELSE
         IF (BC_2L == 0) phi(0:(N1+1),S2p-1,0:(N3+1)) = ghost2UR(0:(N1+1),0:(N3+1))
         IF (BC_2U == 0) phi(0:(N1+1),N2p+1,0:(N3+1)) = ghost2LR(0:(N1+1),0:(N3+1))
      END IF
      
      !-------------------------------------------------------------------------------------------------------------------
      
      IF (ss == 1) THEN
         IF (BC_2L == -1) THEN
            phi(0:(N1+1),S2p  ,0:(N3+1)) = phi(0:(N1+1),S2p,0:(N3+1)) + phi(0:(N1+1),N2p+1,0:(N3+1))
            phi(0:(N1+1),N2p  ,0:(N3+1)) = phi(0:(N1+1),N2p,0:(N3+1)) + phi(0:(N1+1),S2p-1,0:(N3+1))
         END IF
      ELSE
         IF (BC_2L == -1) THEN
            phi(0:(N1+1),S2p-1,0:(N3+1)) = phi(0:(N1+1),N2p,0:(N3+1))
            phi(0:(N1+1),N2p+1,0:(N3+1)) = phi(0:(N1+1),S2p,0:(N3+1))
         END IF
      END IF
      
      !-------------------------------------------------------------------------------------------------------------------
      
      IF (vel_dir == dir) THEN
         IF (BC_2L .GT. 0) phi(0:(N1+1),  -1,0:(N3+1)) = 0.
         IF (BC_2U .GT. 0) phi(0:(N1+1),N2+1,0:(N3+1)) = 0.
         
         IF (BC_2L  == -2) phi(0:(N1+1),   0,0:(N3+1)) = -phi(0:(N1+1),   1,0:(N3+1))
         IF (BC_2U  == -2) phi(0:(N1+1),N2  ,0:(N3+1)) = -phi(0:(N1+1),N2-1,0:(N3+1))
      ELSE
         IF (BC_2L .GT. 0) phi(0:(N1+1),   0,0:(N3+1)) = 0.
         IF (BC_2U .GT. 0) phi(0:(N1+1),N2+1,0:(N3+1)) = 0.
         
         !IF (ss == 1) THEN
         !   IF (BC_2L == -2) phi(0:(N1+1),1 ,0:(N3+1)) = 2.*phi(0:(N1+1),1 ,0:(N3+1))
         !   IF (BC_2U == -2) phi(0:(N1+1),N2,0:(N3+1)) = 2.*phi(0:(N1+1),N2,0:(N3+1))
         !ELSE
            IF (BC_2L == -2) phi(0:(N1+1),   0,0:(N3+1)) = 0.
            IF (BC_2U == -2) phi(0:(N1+1),N2+1,0:(N3+1)) = 0.
         !END IF
      END IF
      
   END IF
   
   !======================================================================================================================
   
   IF (dir == 3 .AND. dimens == 3) THEN
      
      length = (N1+2)*(N2+2)
      
      IF (BC_3L == 0) CALL MPI_IRECV(ghost3UR,length,MPI_REAL8,rank3L,5,COMM_CART,req3L,merror)
      IF (BC_3U == 0) CALL MPI_IRECV(ghost3LR,length,MPI_REAL8,rank3U,6,COMM_CART,req3U,merror)
      
      IF (BC_3U == 0) ghost3US(0:(N1+1),0:(N2+1)) = phi(0:(N1+1),0:(N2+1),N3p+ss)
      IF (BC_3L == 0) ghost3LS(0:(N1+1),0:(N2+1)) = phi(0:(N1+1),0:(N2+1),S3p-ss)
      
      IF (BC_3U == 0) CALL MPI_SEND(ghost3US,length,MPI_REAL8,rank3U,5,COMM_CART,merror)
      IF (BC_3L == 0) CALL MPI_SEND(ghost3LS,length,MPI_REAL8,rank3L,6,COMM_CART,merror)
      
      IF (BC_3L == 0) CALL MPI_WAIT(req3L,status,merror)
      IF (BC_3U == 0) CALL MPI_WAIT(req3U,status,merror)
      
      IF (BC_3L == 0) CALL pseudocall(ghost3UR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
      IF (BC_3U == 0) CALL pseudocall(ghost3LR)
      
      IF (ss == 1) THEN
         IF (BC_3L == 0) phi(0:(N1+1),0:(N2+1),S3p  ) = ghost3UR(0:(N1+1),0:(N2+1)) + phi(0:(N1+1),0:(N2+1),S3p)
         IF (BC_3U == 0) phi(0:(N1+1),0:(N2+1),N3p  ) = ghost3LR(0:(N1+1),0:(N2+1)) + phi(0:(N1+1),0:(N2+1),N3p)
      ELSE
         IF (BC_3L == 0) phi(0:(N1+1),0:(N2+1),S3p-1) = ghost3UR(0:(N1+1),0:(N2+1))
         IF (BC_3U == 0) phi(0:(N1+1),0:(N2+1),N3p+1) = ghost3LR(0:(N1+1),0:(N2+1))
      END IF
      
      !-------------------------------------------------------------------------------------------------------------------
      
      IF (ss == 1) THEN
         IF (BC_3L == -1) THEN
            phi(0:(N1+1),0:(N2+1),S3p  ) = phi(0:(N1+1),0:(N2+1),S3p) + phi(0:(N1+1),0:(N2+1),N3p+1)
            phi(0:(N1+1),0:(N2+1),N3p  ) = phi(0:(N1+1),0:(N2+1),N3p) + phi(0:(N1+1),0:(N2+1),S3p-1)
         END IF
      ELSE
         IF (BC_3L == -1) THEN
            phi(0:(N1+1),0:(N2+1),S3p-1) = phi(0:(N1+1),0:(N2+1),N3p)
            phi(0:(N1+1),0:(N2+1),N3p+1) = phi(0:(N1+1),0:(N2+1),S3p)
         END IF
      END IF
      
      !-------------------------------------------------------------------------------------------------------------------
      
      IF (vel_dir == dir) THEN
         IF (BC_3L .GT. 0) phi(0:(N1+1),0:(N2+1),  -1) = 0.
         IF (BC_3U .GT. 0) phi(0:(N1+1),0:(N2+1),N3+1) = 0.
         
         IF (BC_3L  == -2) phi(0:(N1+1),0:(N2+1),   0) = -phi(0:(N1+1),0:(N2+1),   1)
         IF (BC_3U  == -2) phi(0:(N1+1),0:(N2+1),N3  ) = -phi(0:(N1+1),0:(N2+1),N3-1)
      ELSE
         IF (BC_3L .GT. 0) phi(0:(N1+1),0:(N2+1),   0) = 0.
         IF (BC_3U .GT. 0) phi(0:(N1+1),0:(N2+1),N3+1) = 0.
         
         !IF (ss == 1) THEN
         !   IF (BC_3L == -2) phi(0:(N1+1),0:(N2+1),1 ) = 2.*phi(0:(N1+1),0:(N2+1),1 )
         !   IF (BC_3U == -2) phi(0:(N1+1),0:(N2+1),N3) = 2.*phi(0:(N1+1),0:(N2+1),N3)
         !ELSE
            IF (BC_3L == -2) phi(0:(N1+1),0:(N2+1),   0) = 0.
            IF (BC_3U == -2) phi(0:(N1+1),0:(N2+1),N3+1) = 0.
         !END IF
      END IF
      
   END IF
   
   !======================================================================================================================
   
   
   END SUBROUTINE exchange_part
  
  
  
  
END MODULE mod_exchange
