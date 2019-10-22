#ifndef MOD_SETUP_WRAPPER_H
#define MOD_SETUP_WRAPPER_H
// function prototype/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_setup_MOD_init_general(void);
	void __mod_setup_MOD_init_parallel(void);
	void __mod_setup_MOD_init_boundaries(void);
	void __mod_setup_MOD_init_limits(void);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_init_general(){__mod_setup_MOD_init_general();}
inline void mpt_init_parallel(){__mod_setup_MOD_init_parallel();}
inline void mpt_init_boundaries(){__mod_setup_MOD_init_boundaries();}
inline void mpt_init_limits(){__mod_setup_MOD_init_limits();}
#endif //MOD_SETUP_WRAPPER_H
