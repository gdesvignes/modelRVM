#pragma once
extern "C" void polychord_c_interface(
				      double (*)(double*,int,double*,int,void*), 
        void (*)(double*,double*,int,void*), 
        void (*)(int,int,int,double*,double*,double*,double,double), 
        int,
        int,
        int,
        int,
        bool,
        int,
        double,
        double,
        int,
        double,
        bool,
        bool,
        bool,
        bool,
        bool,
        bool,
        bool,
        bool,
        bool,
        bool,
        bool,
        double,
        int,
        int,
        char*,
        char*,
        int,
        double*,
        int*,
        int,
        double*,
        int*,
        int,
#ifdef USE_MPI
		  MPI_Fint&
#else
		  int&
#endif
		  );

extern "C" void polychord_c_interface_ini(
					  double (*)(double*,int,double*,int,void*), 
        void (*)(), 
        char*,
#ifdef USE_MPI
		  MPI_Fint&
#else
		  int&
#endif
        );
