#ifndef MOD_COEFFS_WRAPPER_H
#define MOD_COEFFS_WRAPPER_H
#include<cstring>
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_coeffs_MOD_fd_coeffs(void);
	void __mod_coeffs_MOD_test_diff(int* dir, int* abl,
			double* xc, double* xe, int* nmax, int* bl, int* bu,
			double* cc, int* bcl, int* grid_type, const char *name,
			int name_length);
	void __mod_coeffs_MOD_diff_coeffs(int* dir, int* abl, int* upwind, int* mapping_yes,
			int* dim_ncb, int* n_coeff_bound,
			double* xc, double* xe, int* bcl, int* bcu, int* grid_type,
			int* nmax, int* bl, int* bu, double* cc);
	void __mod_coeffs_MOD_diff_coeffs_exact(int* abl, int* n_coeff, double* deltax,
			double* cc);
	void __mod_coeffs_MOD_test_coeffs(void);
	void __mod_coeffs_MOD_fd_coeffs_solver(int* abl, int* filter_yes, int* n_coeff,
			double* deltax, double* cc);
	void __mod_coeffs_MOD_fd_coeffs_solver_integral(int* n_coeff,
			double* deltax, double* dxl, double* dxu, double* cc);
	void __mod_coeffs_MOD_interp_coeffs(void);
	void __mod_coeffs_MOD_interp_coeffs_helm(void);
	void __mod_coeffs_MOD_restr_coeffs(void);
	void __mod_coeffs_MOD_restr_coeffs_helm(void);
	void __mod_coeffs_MOD_get_stencil(void);
	void __mod_coeffs_MOD_get_stencil_helm(void);
	void __mod_coeffs_MOD_get_stencil_transp(void);
	void __mod_coeffs_MOD_get_weights(void);
	void __mod_coeffs_MOD_matrix_invert(int* n, double* matrix, double* matrix_inv);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_fd_coeffs(){__mod_coeffs_MOD_fd_coeffs();}
inline void mpt_test_diff(int& dir, int& abl, 
		double* xc, double* xe, int& nmax, int& bl, int& bu, 
		double* cc, int& bcl, int& grid_type, const char* name){
	__mod_coeffs_MOD_test_diff(&dir, &abl, xc, xe, &nmax, &bl, &bu, cc, &bcl, 
			&grid_type, name, strlen(name));
}
inline void mpt_diff_coeffs(int& dir, int& abl, int& upwind, int& mapping_yes,
		int& dim_ncb, int* n_coeff_bound, 
		double* xc, double* xe, int& bcl, int& bcu, int& grid_type,
		int& nmax, int& bl, int& bu, double* cc){
	__mod_coeffs_MOD_diff_coeffs(&dir, &abl, &upwind, &mapping_yes,
			&dim_ncb, n_coeff_bound, xc, xe, &bcl, &bcu, &grid_type,
			&nmax, &bl, &bu, cc);
}
inline void mpt_diff_coeffs_exact(int& abl, int& n_coeff, double& deltax, double* cc){
	__mod_coeffs_MOD_diff_coeffs_exact(&abl, &n_coeff, &deltax, cc);
}
inline void mpt_test_coeffs(){__mod_coeffs_MOD_test_coeffs();}
inline void mpt_fd_coeffs_solver(int& abl, int& filter_yes, int& n_coeff, double* deltax,
		double* cc){
	__mod_coeffs_MOD_fd_coeffs_solver(&abl, &filter_yes, &n_coeff, deltax, cc);
}
inline void mpt_fd_coeffs_solver_integral(int& n_coeff, double* deltax, double& dxl,
		double& dxu, double* cc){
	__mod_coeffs_MOD_fd_coeffs_solver_integral(&n_coeff, deltax, &dxl, &dxu, cc);
}
inline void mpt_interp_coeffs(){__mod_coeffs_MOD_interp_coeffs();}
inline void mpt_interp_coeffs_helm(){__mod_coeffs_MOD_interp_coeffs_helm();}
inline void mpt_restr_coeffs(){__mod_coeffs_MOD_restr_coeffs();}
inline void mpt_restr_coeffs_helm(){__mod_coeffs_MOD_restr_coeffs_helm();}
inline void mpt_get_stencil(){__mod_coeffs_MOD_get_stencil();}
inline void mpt_get_stencil_helm(){__mod_coeffs_MOD_get_stencil_helm();}
inline void mpt_get_stencil_transp(){__mod_coeffs_MOD_get_stencil_transp();}
inline void mpt_get_weights(){__mod_coeffs_MOD_get_weights();}
inline void mpt_matrix_invert(int& n, double* matrix, double* matrix_inv){
	__mod_coeffs_MOD_matrix_invert(&n, matrix, matrix_inv);
}
#endif //MOD_COEFFS_WRAPPER_H
