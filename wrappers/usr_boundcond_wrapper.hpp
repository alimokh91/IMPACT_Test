#ifndef USR_BOUNDCOND_WRAPPER_H
#define USR_BOUNDCOND_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void boundary_vel_stat_(void);
       void boundary_vel_tint_(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_boundary_vel_stat(){boundary_vel_stat_();}
inline void mpt_boundary_vel_tint(){boundary_vel_tint_();}

#endif //USR_BOUNDCOND_WRAPPER_H
