#ifndef USR_INITCOND_WRAPPER_H
#define USR_INITCOND_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void initial_conditions_vel_(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_initial_conditions_vel(){initial_conditions_vel_();}

#endif //USR_INITCOND_WRAPPER_H
