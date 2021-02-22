#ifndef MOD_HELMHOLTZ_WRAPPER_H
#define MOD_HELMHOLTZ_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_helmholtz_MOD_product_helmholtz(double* phi, double* hel);
	void __mod_helmholtz_MOD_product_helmholtz_relax(double* phi, double* hel);
	void __mod_helmholtz_MOD_product_helmholtz_relax_coarse(int* g, double* phi, double* hel);
	void __mod_helmholtz_MOD_relaxation_helmholtz(int* init_yes, int* n_relax, 
			double* bb, double* rel);
	void __mod_helmholtz_MOD_relaxation_helmholtz_coarse(int* init_yes, int* n_relax, int* g,
			double* bb, double* rel);
#ifdef __cplusplus
}
#endif

// function prototypes
inline void mpt_product_helmholtz(double* phi, double* hel){
	__mod_helmholtz_MOD_product_helmholtz(phi, hel);
}
inline void mpt_product_helmholtz_relax(double* phi, double* hel){
	__mod_helmholtz_MOD_product_helmholtz_relax(phi, hel);
}
inline void mpt_product_helmholtz_relax_coarse(int& g, double* phi, double* hel){
	__mod_helmholtz_MOD_product_helmholtz_relax_coarse(&g, phi, hel);
}
inline void mpt_relaxation_helmholtz(int& init_yes, int& n_relax, double* bb, double* rel){
	__mod_helmholtz_MOD_relaxation_helmholtz(&init_yes, &n_relax, bb, rel);
}
inline void mpt_relaxation_helmholtz_coarse(int& init_yes, int& n_relax, int& g, 
		double* bb, double* rel){
	__mod_helmholtz_MOD_relaxation_helmholtz_coarse(&init_yes, &n_relax, &g, bb, rel);
}
#endif //MOD_HELMHOLTZ_WRAPPER_H
