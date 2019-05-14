#ifndef USR_GEOMETRY_WRAPPER_H
#define USR_GEOMETRY_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void get_coords_(void);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void get_coords(){get_coords_();}
#endif //USR_GEOMETRY_WRAPPER_H
