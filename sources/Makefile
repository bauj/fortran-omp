myExe : main.f90
	gfortran -o myExe.exe main.f90 -fopenmp -fno-range-check -I/usr/include/ \
		-L/lib64/ -lfftw3_threads -lfftw3
clean:
	rm -rf myExe.exe
