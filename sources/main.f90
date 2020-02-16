PROGRAM MainProgramTest

	CALL PARALLEL_HELLO_WORLD()
	CALL PARALLEL_COUNT() 

END PROGRAM MainProgramTest      

!____________________________________________________________________________
SUBROUTINE FFTW_REAL_2_COMPLEX
USE, INTRINSIC :: iso_c_binding
INCLUDE 'fftw3.f03'

	RETURN

END SUBROUTINE FFTW_REAL_2_COMPLEX


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
