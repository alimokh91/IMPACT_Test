#ifndef MOD_TIMEINT_WRAPPER_H
#define MOD_TIMEINT_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void __mod_timeint_MOD_timeintegration(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_timeintegration(){__mod_timeint_MOD_timeintegration();}

#endif //MOD_TIMEINT_WRAPPER_H
