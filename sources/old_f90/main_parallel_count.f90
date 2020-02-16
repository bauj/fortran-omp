PROGRAM Parallel_Hello_World
USE OMP_LIB

INTEGER(kind=8) :: partial_Sum, total_Sum,i
REAL(kind=16)   :: o_cpu_time1, o_cpu_time2    

o_cpu_time1 = omp_get_wtime()

!$OMP PARALLEL PRIVATE(partial_Sum) SHARED(total_Sum)
    partial_Sum = 0;
    total_Sum = 0;

    !$OMP DO
    DO i=1,2000000000
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

END PROGRAM Parallel_Hello_World
