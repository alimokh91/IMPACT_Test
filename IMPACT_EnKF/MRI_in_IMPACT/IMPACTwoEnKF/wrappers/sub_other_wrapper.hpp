#ifndef SUB_OTHER_WRAPPER_H
#define SUB_OTHER_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void pseudocall_(double* phi);
	void pseudocall_int_(int* phi);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_pseudocall(double* phi){pseudocall_(phi);}
inline void mpt_pseudocall_int(int* phi){pseudocall_int_(phi);}
#endif //SUB_OTHER_WRAPPER_H
