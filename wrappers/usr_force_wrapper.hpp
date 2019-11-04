#ifndef USR_FORCE_WRAPPER_H
#define USR_FORCE_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void forcing_vel_(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_forcing_vel(){forcing_vel_();}

#endif //USR_FORCE_WRAPPER_H
