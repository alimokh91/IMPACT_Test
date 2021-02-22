#ifndef MOD_LIB_WRAPPER_H
#define MOD_LIB_WRAPPER_H
#include<cstring>
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_lib_MOD_get_dtime(void);
	void __mod_lib_MOD_get_beta(void);
	void __mod_lib_MOD_init_bc(void);
	void __mod_lib_MOD_fill_corners(double* phi);
	void __mod_lib_MOD_level_pressure(void);
	void __mod_lib_MOD_init_alarm(void);
	void __mod_lib_MOD_check_alarm(void);
	void __mod_lib_MOD_check_signal(void);
	void __mod_lib_MOD_num_to_string(int* n_digits, int* num, char* num_char,
			int num_char_length);
	void __mod_lib_MOD_iteration_stats(void);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_get_dtime(){__mod_lib_MOD_get_dtime();}
inline void mpt_get_beta(){__mod_lib_MOD_get_beta();}
inline void mpt_init_bc(){__mod_lib_MOD_init_bc();}
inline void mpt_fill_corners(double* phi){__mod_lib_MOD_fill_corners(phi);}
inline void mpt_level_pressure(){__mod_lib_MOD_level_pressure();}
inline void mpt_init_alarm(){__mod_lib_MOD_init_alarm();}
inline void mpt_check_alarm(){__mod_lib_MOD_check_alarm();}
inline void mpt_check_signal(){__mod_lib_MOD_check_signal();}
inline void mpt_num_to_string(int& n_digits, int& num, char* num_char){
	__mod_lib_MOD_num_to_string(&n_digits, &num, num_char, strlen(num_char));
}
inline void mpt_iteration_stats(){__mod_lib_MOD_iteration_stats();}
#endif //MOD_LIB_WRAPPER_H
