#ifndef USR_CONFIG_WRAPPER_H
#define USR_CONFIG_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void configuration_(void);
#ifdef __cplusplus
}
#endif

// wrapper functions
inline void mpt_configuration(){configuration_();}
#endif //USR_CONFIG_WRAPPER_H
