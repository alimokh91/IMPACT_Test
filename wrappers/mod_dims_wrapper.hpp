#ifndef MOD_DIMS_WRAPPER_H
#define MOD_DIMS_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	// fortran module variables
	// domain and block specifications
	// block numbers
	extern int __mod_dims_MOD_nb1;
	extern int __mod_dims_MOD_nb2;
	extern int __mod_dims_MOD_nb3;
	// entire domain
	extern int __mod_dims_MOD_m1;
	extern int __mod_dims_MOD_m2;
	extern int __mod_dims_MOD_m3;
#ifdef __cplusplus
}
#endif

// variable macros
#define mpt_nb1 __mod_dims_MOD_nb1
#define mpt_nb2 __mod_dims_MOD_nb2
#define mpt_nb3 __mod_dims_MOD_nb3
#define mpt_m1 __mod_dims_MOD_m1
#define mpt_m2 __mod_dims_MOD_m2
#define mpt_m3 __mod_dims_MOD_m3

#endif //MOD_DIMS_WRAPPER_H
