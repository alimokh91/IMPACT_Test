#ifndef USR_POST_WRAPPER_H
#define USR_POST_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void postprocess_(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_postprocess(){postprocess_();}

#endif //USR_POST_WRAPPER_H
