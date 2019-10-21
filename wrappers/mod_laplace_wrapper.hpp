#ifndef MOD_LAPLACE_WRAPPER_H
#define MOD_LAPLACE_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_laplace_MOD_product_div_grad(double* phi, double* lap);
	void __mod_laplace_MOD_product_div_grad_transp(double* phi, double* lap);
	void __mod_laplace_MOD_product_div_grad_relax(int* g, double* phi, double* lap);
	void __mod_laplace_MOD_relaxation_div_grad(int* init_yes, int* n_relax, int* g,
			double* bb, double* rel);
	void __mod_laplace_MOD_relaxation_div_grad_inv(int* init_yes, int* n_relax, int* g,
			double* bb, double* rel);
	void __mod_laplace_MOD_handle_corner_lap(int* g, double* phi);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_product_div_grad(double* phi, double* lap){
	__mod_laplace_MOD_product_div_grad(phi, lap);
}
inline void mpt_product_div_grad_transp(double* phi, double* lap){
	__mod_laplace_MOD_product_div_grad_transp(phi, lap);
}
inline void mpt_product_div_grad_relax(int& g, double* phi, double* lap){
	__mod_laplace_MOD_product_div_grad_relax(&g, phi, lap);
}
inline void mpt_relaxation_div_grad(int& init_yes, int& n_relax, int& g,
		double* bb, double* rel){
	__mod_laplace_MOD_relaxation_div_grad(&init_yes, &n_relax, &g, bb, rel);
}
inline void mpt_relaxation_div_grad_inv(int& init_yes, int& n_relax, int&g, 
		double* bb, double* rel){
	__mod_laplace_MOD_relaxation_div_grad_inv(&init_yes, &n_relax, &g, bb, rel);
}
inline void mpt_handle_corner_lap(int& g, double* phi){
	__mod_laplace_MOD_handle_corner_lap(&g, phi);
}
#endif //MOD_LAPLACE_WRAPPER_H
