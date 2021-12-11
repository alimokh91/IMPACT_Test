#ifndef MOD_TEST_WRAPPER_H
#define MOD_TEST_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif

//fortran subroutines 

       void __mod_test_MOD_test_parameter(void);
       void __mod_test_MOD_test_divergence(void);
       void __mod_test_MOD_test_momentum(void);
       void __mod_test_MOD_analyze_matrix(int* gridtype);


#ifdef __cplusplus
}
#endif

// function wrappers 

inline void mpt_test_parameter(){__mod_test_MOD_test_parameter();}
inline void mpt_test_divergence(){__mod_test_MOD_test_divergence();}
inline void mpt_test_momentum(){__mod_test_MOD_test_momentum();}
inline void mpt_analyze_matrix(int& gridtype){__mod_test_MOD_analyze_matrix(&gridtype);}

#endif //MOD_TEST_WRAPPER_H
