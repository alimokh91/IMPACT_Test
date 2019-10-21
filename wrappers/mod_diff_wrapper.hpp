#ifndef MOD_DIFF_WRAPPER_H
#define MOD_DIFF_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_diff_MOD_divergence(int* m, double* phi, double* div);
	void __mod_diff_MOD_divergence2(double* phi, double* div);
	void __mod_diff_MOD_divergence_transp(int* m, double* phi, double* div);
	void __mod_diff_MOD_gradient(int* m, double* phi, double* grad);
	void __mod_diff_MOD_gradient_transp(int* m, double* phi, double* grad);
	void __mod_diff_MOD_helmholtz(int* m, int* exch_yes, double* phi, double* lap);
	void __mod_diff_MOD_helmholtz_explicit(int* exch_yes);
	void __mod_diff_MOD_nonlinear(int* exch_yes);
	void __mod_diff_MOD_interpolate_vel(int* exch_yes);
	void __mod_diff_MOD_outflow_bc(void);
	void __mod_diff_MOD_bc_extrapolation(int* m, double* phi);
	void __mod_diff_MOD_bc_extrapolation_transp(int* m, double* phi);
	void __mod_diff_MOD_interpolate_pre_vel(int* exch_yes, int* m, 
			int* ss1, int* ss2, int* ss3, 
			int* nn1, int* nn2, int* nn3, 
			double* phi, double* inter);
	void __mod_diff_MOD_interpolate_vel_pre(int* exch_yes, int* m, 
			int* ss1, int* ss2, int* ss3, 
			int* nn1, int* nn2, int* nn3, 
			double* phi, double* inter);
	void __mod_diff_MOD_interpolate2_pre_vel(int* exch_yes, int* m, 
			double* phi, double* inter);
	void __mod_diff_MOD_interpolate2_vel_pre(int* exch_yes, int* m, 
			double* phi, double* inter);
	void __mod_diff_MOD_first_adv_pre(int* exch_yes, int* m,
			int* ss1, int* ss2, int* ss3, 
			int* nn1, int* nn2, int* nn3, 
			double* phi, double* der, double* adv,
			int* upwind_yes);
	void __mod_diff_MOD_first_adv_vel(int* exch_yes, int* m,
			int* ss1, int* ss2, int* ss3, 
			int* nn1, int* nn2, int* nn3, 
			double* phi, double* der, double* adv,
			int* upwind_yes);
	void __mod_diff_MOD_helmholtz_pre_explicit(int* exch_yes, double* phi, double* lap);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_divergence(int& m, double* phi, double* div){
	__mod_diff_MOD_divergence(&m, phi, div);
}
inline void mpt_divergence2(double* phi, double* div){
	__mod_diff_MOD_divergence2(phi, div);
}
inline void mpt_divergence_transp(int& m, double* phi, double* div){
	__mod_diff_MOD_divergence_transp(&m, phi, div);
}
inline void mpt_gradient(int& m, double* phi, double* grad){
	__mod_diff_MOD_gradient(&m, phi, grad);
}
inline void mpt_gradient_transp(int& m, double* phi, double* grad){
	__mod_diff_MOD_gradient_transp(&m, phi, grad);
}
inline void mpt_helmholtz(int& m, int& exch_yes, double* phi, double* lap){
	__mod_diff_MOD_helmholtz(&m, &exch_yes, phi, lap);
}
inline void mpt_helmholtz_explicit(int& exch_yes){
	__mod_diff_MOD_helmholtz_explicit(&exch_yes);
}
inline void mpt_nonlinear(int& exch_yes){
	__mod_diff_MOD_nonlinear(&exch_yes);
}
inline void mpt_interpolate_vel(int exch_yes){
	__mod_diff_MOD_interpolate_vel(&exch_yes);
}
inline void mpt_outflow_bc(){
	__mod_diff_MOD_outflow_bc();
}
inline void mpt_bc_extrapolation(int& m, double* phi){
	__mod_diff_MOD_bc_extrapolation(&m, phi);
}
inline void mpt_bc_extrapolation_transp(int& m, double* phi){
	__mod_diff_MOD_bc_extrapolation_transp(&m, phi);
}
inline void mpt_interpolate_pre_vel(int& exch_yes, int& m, 
		int& ss1, int& ss2, int& ss3,
		int& nn1, int& nn2, int& nn3,
		double* phi, double* inter){
	__mod_diff_MOD_interpolate_pre_vel(&exch_yes, &m,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, phi, inter);
}
inline void mpt_interpolate_vel_pre(int& exch_yes, int& m, 
		int& ss1, int& ss2, int& ss3,
		int& nn1, int& nn2, int& nn3,
		double* phi, double* inter){
	__mod_diff_MOD_interpolate_vel_pre(&exch_yes, &m,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, phi, inter);
}
inline void mpt_interpolate2_pre_vel(int& exch_yes, int& m,
		double* phi, double* inter){
	__mod_diff_MOD_interpolate2_pre_vel(&exch_yes, &m, phi, inter);
}
inline void mpt_interpolate2_vel_pre(int& exch_yes, int& m,
		double* phi, double* inter){
	__mod_diff_MOD_interpolate2_vel_pre(&exch_yes, &m, phi, inter);
}
inline void mpt_first_adv_pre(int& exch_yes, int& m, 
		int& ss1, int& ss2, int& ss3,
		int& nn1, int& nn2, int& nn3,
		double* phi, double* der, double* adv,
		int& upwind_yes){
	__mod_diff_MOD_first_adv_pre(&exch_yes, &m,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, phi, der, adv,
			&upwind_yes);
}
inline void mpt_first_adv_vel(int& exch_yes, int& m, 
		int& ss1, int& ss2, int& ss3,
		int& nn1, int& nn2, int& nn3,
		double* phi, double* der, double* adv,
		int& upwind_yes){
	__mod_diff_MOD_first_adv_vel(&exch_yes, &m,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, phi, der, adv,
			&upwind_yes);
}
inline void mpt_helmholtz_pre_explicit(int& exch_yes, double* phi, double* lap){
	__mod_diff_MOD_helmholtz_pre_explicit(&exch_yes, phi, lap);
}
#endif //MOD_DIFF_WRAPPER_H
