PROGRAM ordered_hello_world_omp
USE OMP_LIB

INTEGER(kind=8) :: thread_id, nb_threads,i 

!$OMP PARALLEL PRIVATE(thread_id)

	nb_threads = OMP_GET_MAX_THREADS()
	thread_id  = OMP_GET_THREAD_NUM()

	DO i = 0, nb_threads-1
		IF (i == thread_id) THEN
			PRINT*, "Hello world from process :", thread_id
		END IF
		!$OMP BARRIER
	END DO

!$OMP END PARALLEL

END PROGRAM ordered_hello_world_omp


