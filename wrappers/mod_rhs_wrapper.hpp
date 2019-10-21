#ifndef MOD_RHS_WRAPPER_H
#define MOD_RHS_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void __mod_rhs_MOD_rhs_vel(void);
#ifdef __cplusplus
}
#endif
// function wrappers

inline void mpt_rhs_vel(){__mod_rhs_MOD_rhs_vel();}

#endif //MOD_RHS_WRAPPRE_H

