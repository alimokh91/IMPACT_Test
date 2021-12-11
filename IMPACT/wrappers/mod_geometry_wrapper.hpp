#ifndef MOD_GEOMETRY_WRAPPER_H
#define MOD_GEOMETRY_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif

       void __mod_geometry_MOD_coordinates(void);

#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_coordinates(){__mod_geometry_MOD_coordinates();}

#endif //MOD_GEOMETRY_WRAPPER_H
