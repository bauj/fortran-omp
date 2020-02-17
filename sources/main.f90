PROGRAM MainProgramTest
	
	USE OMP_LIB
	IMPLICIT NONE

	REAL(kind=8),    ALLOCATABLE :: real_init(:,:,:), real_recons(:,:,:)
	COMPLEX(kind=8), ALLOCATABLE :: tf_init(:,:,:) 
	INTEGER(kind=4) :: i,j,k,nx,ny,nz,inputs(4)
	REAL(kind=8)    :: pi, pi2, val1, val2, val3
	REAL(kind=8)    :: otime1, otime2, time1, time2 

	inputs = READ_INPUT("input.data")

	!$ CALL OMP_SET_NUM_THREADS(inputs(1))

	CALL cpu_time(time1)
	!$ otime1 = OMP_GET_WTIME()
	! Parameters
	nx = inputs(2)
	ny = inputs(3)
	nz = inputs(4)

	pi  = acos(-1.d0)
	pi2 = 2.d0 * pi

	ALLOCATE(real_init(nx,ny,nz))
	ALLOCATE(real_recons(nx,ny,nz))
	ALLOCATE(tf_init(nx/2+1,ny,nz))
  	
	real_init(:,:,:)   = 0.d0
	real_recons(:,:,:) = 0.d0
	tf_init(:,:,:)     = 0.d0

	!$OMP PARALLEL PRIVATE(k,j,i) SHARED(real_init,tf_init,real_recons,nx,ny,nz) 

	!$OMP DO 	
	DO k=1,nz
	DO j=1,ny
	DO i=1,nx
		real_init(i,j,k) = COS(i*pi2/nx)
	END DO
	END DO
	END DO
	!$OMP END DO

	PRINT*, "Entering fftw ..."
	CALL FFTW_REAL_2_COMPLEX(nx,ny,nz,real_init,tf_init)

	CALL FFTW_COMPLEX_2_REAL(nx,ny,nz,tf_init,real_recons)

    !$OMP END PARALLEL
    
    val1 = MAXVAL(ABS(tf_init))
    val2 = MAXVAL(ABS(real_init(:,:,:)-real_recons(:,:,:)))
    PRINT*, "Maxval modulo tf init =", val1
    PRINT*, "Maxval erreur reconst =", val2
    DEALLOCATE(real_init)
    DEALLOCATE(real_recons)
    DEALLOCATE(tf_init)
    
    CALL cpu_time(time2)
    !$ otime2 = OMP_GET_WTIME()
    
    PRINT*, "Elapsed time : ",otime2-otime1 ,"s."
    
!	CALL PARALLEL_HELLO_WORLD()
!	CALL PARALLEL_COUNT() 

CONTAINS

	FUNCTION READ_INPUT(filenm) RESULT(inputs)

		IMPLICIT NONE

		CHARACTER(10), INTENT(IN)   :: filenm
		INTEGER(kind=4)             :: inputs(4)
		CHARACTER(LEN=8), PARAMETER :: FMT1 = '(a16,i3)'
	    CHARACTER(16)               :: inputStr
        
	    OPEN(UNIT=2, FILE=filenm)
	    READ(2,FMT1) inputStr, inputs(1)
	    READ(2,FMT1) inputStr, inputs(2)
	    READ(2,FMT1) inputStr, inputs(3)
	    READ(2,FMT1) inputStr, inputs(4)
        CLOSE(UNIT=2)
		
		RETURN 

	END FUNCTION READ_INPUT

END PROGRAM MainProgramTest      

!____________________________________________________________________________
SUBROUTINE FFTW_REAL_2_COMPLEX(lx,ly,lz,tab_in, tab_tf)

    USE OMP_LIB
	IMPLICIT NONE

	INTEGER(kind=4), INTENT(IN)  :: lx,ly,lz	
	REAL(kind=8),    INTENT(IN)  :: tab_in(lx,ly,lz)
	COMPLEX(kind=8), INTENT(OUT) :: tab_tf(lx/2+1, ly, lz)
	INTEGER(kind=8)	:: plan
	INTEGER(kind=8)	:: N_threads 
	REAL(kind=8)	:: iRet 
	INTEGER				:: fftw_estimate
	PARAMETER(fftw_estimate=64)


	!$OMP SINGLE
    N_threads = OMP_GET_MAX_THREADS()
	
	PRINT*, "Number of threads : ", N_threads	
	CALL dfftw_init_threads(iRet)
	IF (iRet==0) THEN
		PRINT*, "Problem initializing FFTW mthreaded ... Quitting ..."
		STOP
	END IF

    CALL dfftw_plan_with_nthreads(N_threads)

	CALL dfftw_plan_dft_r2c_3d(plan,lx,ly,lz,tab_in,tab_tf,fftw_estimate)

	CALL dfftw_execute_dft_r2c(plan,tab_in,tab_tf)

	tab_tf(:,:,:) = tab_tf(:,:,:)/(REAL(lx)*REAL(ly)*REAL(lz))

	CALL dfftw_destroy_plan(plan)

	CALL dfftw_cleanup_threads()

	!$OMP END SINGLE
	RETURN

END SUBROUTINE FFTW_REAL_2_COMPLEX

!____________________________________________________________________________
SUBROUTINE FFTW_COMPLEX_2_REAL(lx,ly,lz,tab_tf, tab_out)

    USE OMP_LIB
	IMPLICIT NONE
    
    INTEGER(kind=4), INTENT(IN)  :: lx,ly,lz	
    COMPLEX(kind=8), INTENT(IN)  :: tab_tf(lx/2+1, ly, lz)
    REAL(kind=8),    INTENT(OUT) :: tab_out(lx,ly,lz)
    COMPLEX(kind=8) :: temp_tf(lx/2+1, ly, lz)
    INTEGER(kind=8)	:: plan
    INTEGER(kind=8)	:: N_threads 
    REAL(kind=8)	:: iRet 
    INTEGER			:: fftw_estimate
    PARAMETER(fftw_estimate=64)
    
    !$OMP SINGLE
    temp_tf(:,:,:) = tab_tf(:,:,:)
    N_threads = OMP_GET_MAX_THREADS()
    
    PRINT*, "Number of threads : ", N_threads	
    CALL dfftw_init_threads(iRet)
    IF (iRet==0) THEN
    	PRINT*, "Problem initializing FFTW mthreaded ... Quitting ..."
    	STOP
    END IF
    
    CALL dfftw_plan_with_nthreads(N_threads)
    
    CALL dfftw_plan_dft_c2r_3d(plan,lx,ly,lz,temp_tf,tab_out,fftw_estimate)
    
    CALL dfftw_execute_dft_c2r(plan,temp_tf,tab_out)
    
    CALL dfftw_destroy_plan(plan)
    
    CALL dfftw_cleanup_threads()
    
    !$OMP END SINGLE
    RETURN

END SUBROUTINE FFTW_COMPLEX_2_REAL

!____________________________________________________________________________
SUBROUTINE PARALLEL_HELLO_WORLD
USE OMP_LIB
IMPLICIT NONE

	INTEGER(kind=8) :: thread_id, nb_threads,i 
	
	nb_threads = OMP_GET_MAX_THREADS()
	PRINT*, "Nombre de threads =       ", nb_threads

	!$OMP PARALLEL PRIVATE(thread_id) SHARED(nb_threads)
	
		thread_id  = OMP_GET_THREAD_NUM()
	
		DO i = 0, nb_threads-1
			IF (i == thread_id) THEN
				PRINT*, "Hello world from process :", thread_id
			END IF
			!$OMP BARRIER
		END DO
	
	!$OMP END PARALLEL

	RETURN

END SUBROUTINE PARALLEL_HELLO_WORLD

!____________________________________________________________________________
SUBROUTINE PARALLEL_COUNT()
USE OMP_LIB
IMPLICIT NONE

	INTEGER(kind=8) :: partial_Sum, total_Sum,i
	REAL(kind=16)   :: o_cpu_time1, o_cpu_time2    
	
	o_cpu_time1 = omp_get_wtime()
	
	!$OMP PARALLEL PRIVATE(partial_Sum) SHARED(total_Sum)
	    partial_Sum = 0
	    total_Sum = 0
	
	    !$OMP DO
	    DO i=1,200000000
	        partial_Sum = partial_Sum + i
	    END DO
	    !$OMP END DO
	
	    !$OMP CRITICAL
	        total_Sum = total_Sum + partial_Sum
	    !$OMP END CRITICAL
	
	!$OMP END PARALLEL
	PRINT *, "Total Sum: ", total_Sum
	
	o_cpu_time2 = omp_get_wtime()
	PRINT*, "Elapsed time :",o_cpu_time2-o_cpu_time1,"s."

	RETURN

END SUBROUTINE PARALLEL_COUNT

!____________________________________________________________________________
