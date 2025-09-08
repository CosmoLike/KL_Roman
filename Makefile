class:
	cd ../cosmolike_core/class; $(MAKE)

jpl:
	gcc -shared -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib like_fourier.c -o like_fourier.so -fPIC -lgsl -lfftw3 -lgslcblas -std=gnu99 -O3 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass
	gcc -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
	gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/home/teifler/include -L/home/teifler/lib -o like_fourier like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass

ocelote:
	 gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -shared -o like_fourier.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
	 gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -o like_fourier like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
	 gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass

v2:
	 gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -shared -o like_fourier_v2.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
	 gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -o like_fourier like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
	 gcc -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass


header_home := -I/opt/homebrew/opt/gsl/include -L/opt/homebrew/opt/gsl/lib \
			-I/opt/homebrew/opt/fftw/include -L/opt/homebrew/opt/fftw/lib
FLAGS := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
		-lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 -ffast-math -funroll-loops \
		-std=gnu99 -L../cosmolike_core/class -lclass
home:
	gcc -shared -o like_fourier.so -fPIC like_fourier.c $(header_home) $(FLAGS)
	gcc -o like_fourier like_fourier.c $(header_home) $(FLAGS)
	gcc -o ./compute_covariances_fourier compute_covariances_fourier.c $(header_home) $(FLAGS)
compute_cov:
	gcc -o ./compute_covariances_fourier compute_covariances_fourier.c $(header_home) $(FLAGS)