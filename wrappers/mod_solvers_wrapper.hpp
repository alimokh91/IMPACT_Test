#ifndef MOD_SOLVERS_WRAPPER_H
#define MOD_SOLVERS_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void __mod_solvers_MOD_outer_iteration(void);
       void __mod_solvers_MOD_explicit(void);
       void __mod_solvers_MOD_twostep(void) ;      
       void __mod_solvers_MOD_force_massflow(int* impl_yes);
       void __mod_solvers_MOD_apply_nullspace(double* res);    //inout
       void __mod_solvers_MOD_get_nullspace(void); 
       void __mod_solvers_MOD_solve_nullspace(double* eps, int* g, double* psi, int* problem_type);       
       void __mod_solvers_MOD_solve_helmholtz(int* m, double* epsu, int* n_it_max, int* init_yes,
		       double* bb, double* phi, int* quiet_yes1, int* quiet_yes2); 
       void __mod_solvers_MOD_multigridv(int* init_yes, int* gstart, double* bb, double* phi, int* problem_type);
       void __mod_solvers_MOD_multigridf(int* init_yes, int* gstart, double* bb, double* phi, int* problem_type);
       void __mod_solvers_MOD_restrict(int* add_yes, int* g, double* coarse, double* fine1, double* fine2);
       void __mod_solvers_MOD_interpolate(int* add_yes, int* g, double* coarse, double* fine, double* work);
       void __mod_solvers_MOD_restrict_helmholtz(double* coarse, double* fine1, double* fine2);
       void __mod_solvers_MOD_interpolate_helmholtz(double* coarse, double* fine, double* work);
       void __mod_solvers_MOD_bicgstab(double* eps, int* n_it_max, int* init_yes, 
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3, 
		       double* bb, double* phi, int* problem_type, 
		       int* quiet_yes1, int* quiet_yes2, int* preconditioner);
       void __mod_solvers_MOD_richardson(double* eps, int* n_it_max, int* init_yes, 
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3, 
		       double* bb, double* phi, int* problem_type, 
		       int* quiet_yes1, int* quiet_yes2, int* preconditioner);  
       void __mod_solvers_MOD_get_norms(int* ss1, int* ss2, int* ss3, int* nn1, int* nn2, int *nn3, 
		       double* phi, int* problem_type, int* inf_yes, int* two_yes, double* norminf, double* normtwo);
       void __mod_solvers_MOD_product_scalar(int* ss1,int* ss2, int* ss3, int* nn1, int* nn2,int* nn3, 
		       double* phi1, double* phi2, double* scalar);    
       void __mod_solvers_MOD_multadd1(int* ss1, int* ss2, int* ss3, int* nn1, int* nn2, int* nn3,
		       double* mult, double* vec, double* phi);
       void __mod_solvers_MOD_multadd2(int* ss1, int* ss2, int* ss3, int* nn1, int* nn2, int* nn3, 
		       double* mult1, double* mult2, double* vec, double* phi1, double* phi2, int* init_yes);
       void __mod_solvers_MOD_status_iteration(double* eps, double* norm, int* counter, int* n_restarts,
		       int* exit_yes, int* quiet_yes1, int* quiet_yes2);
       void __mod_solvers_MOD_handle_corner_rhs(int* g, double* phi);
       void __mod_solvers_MOD_bicgstab2(double* eps, int* n_it_max, int* init_yes, 
		       int* n1, int* n2, int* n3, int* g,
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3, 
		       double* bb, double* phi, int* problem_type, int *quiet_yes1,
		       int* quiet_yes2, int* preconditioner);      
       void __mod_solvers_MOD_get_norms2(int* g, int* n1, int* n2, int* n3, 
		       int* ss1, int* ss2, int* ss3, int* nn1, int* nn2, int* nn3,
		       double* phi, int* problem_type, int* inf_yes, int* two_yes, 
		       double* norminf, double* normtwo);
       void __mod_solvers_MOD_product_scalar2(int*g, int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3, double* phi1, double* phi2, double* scalar);   
       void __mod_solvers_MOD_apply_nullspace2(int* g, double* psi, double* res, int* problem_type);
       void __mod_solvers_MOD_relax_restrict(int* init_yes, int* nullspace_yes, int* g, 
		       double* psi, double* bb, double* phi, double* work, double* coarse, int* problem_type);
       void __mod_solvers_MOD_interpolate_relax(int* g,double* bb, double* phi, double* work, double* coarse, int* problem_type);
       void __mod_solvers_MOD_plain_restrict(int* nullspace_yes, int* g, 
		       double* psi, double* bb, double* work, double* coarse, int* problem_type);   
       void __mod_solvers_MOD_interpolate_mg(int* g, double* bb, double* phi, double* work, double* coarse, int* problem_type);
       void __mod_solvers_MOD_relax_bottom(int* init_yes, int* nullspace_yes, int* g, 
		       double* psi, double* bb, double* phi, int* problem_type);   
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_outer_iteration(){__mod_solvers_MOD_outer_iteration();}
inline void mpt_explicit(){__mod_solvers_MOD_explicit();}
inline void mpt_twostep(){__mod_solvers_MOD_twostep();}
inline void mpt_force_massflow(int& impl_yes){__mod_solvers_MOD_force_massflow(&impl_yes);}
inline void mpt_apply_nullspace( double* res){__mod_solvers_MOD_apply_nullspace(res);}
inline void mpt_get_nullspace(){__mod_solvers_MOD_get_nullspace();}
inline void mpt_solve_nullspace(double& eps, int& g, double* psi, int& problem_type){
	__mod_solvers_MOD_solve_nullspace(&eps, &g, psi, &problem_type);
}
inline void mpt_solve_helmholtz(int& m, double& epsu, int& n_it_max, int& init_yes,
		double* bb, double* phi, int& quiet_yes1, int& quiet_yes2){
	__mod_solvers_MOD_solve_helmholtz(&m, &epsu, &n_it_max, &init_yes,
			bb, phi, &quiet_yes1, &quiet_yes2);
}
inline void mpt_multigridv(int& init_yes, int& gstart, double* bb, double* phi, int& problem_type){
	__mod_solvers_MOD_multigridv(&init_yes, &gstart, bb, phi, &problem_type);
}
inline void mpt_multigridf(int& init_yes, int& gstart, double* bb, double* phi, int& problem_type){
	__mod_solvers_MOD_multigridf(&init_yes, &gstart, bb, phi, &problem_type);
}
inline void mpt_restrict(int& add_yes, int& g, double* coarse, double* fine1, double* fine2){
	__mod_solvers_MOD_restrict(&add_yes, &g, coarse, fine1, fine2);
}
inline void mpt_interpolate(int& add_yes, int& g, double* coarse, double* fine, double* work){
	__mod_solvers_MOD_interpolate(&add_yes, &g, coarse, fine, work);
}
inline void mpt_restrict_helmholtz(double* coarse, double*  fine1, double* fine2){
	__mod_solvers_MOD_restrict_helmholtz(coarse, fine1, fine2);
}
inline void mpt_interpolate_Helmholtz(double* coarse, double* fine, double* work){
	__mod_solvers_MOD_interpolate_helmholtz(coarse, fine, work);
}
inline void mpt_bicgstab(double& eps, int& n_it_max, int& init_yes, 
		int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3,
		double* bb, double* phi, int& problem_type, int& quiet_yes1,
		int& quiet_yes2, int& preconditioner){
	__mod_solvers_MOD_bicgstab(&eps, &n_it_max, &init_yes,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3,
			bb, phi, &problem_type, &quiet_yes1, &quiet_yes2, &preconditioner);
}      
inline void mpt_richardson(double& eps, int& n_it_max, int& init_yes, 
		int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3, 
		double* bb, double* phi, int& problem_type, int& quiet_yes1,
		int& quiet_yes2, int& preconditioner){
	__mod_solvers_MOD_richardson(&eps, &n_it_max, &init_yes,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, 
			bb, phi, &problem_type, &quiet_yes1, &quiet_yes2, &preconditioner);
}
inline void mpt_get_norms(int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3, 
		double* phi, int& problem_type, int& inf_yes, int& two_yes, double& norminf, double& normtwo){
	__mod_solvers_MOD_get_norms(&ss1, &ss2, &ss3, &nn1, &nn2, &nn3,
			phi, &problem_type, &inf_yes, &two_yes, &norminf, &normtwo);
}
inline void mpt_product_scalar(int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int nn3, 
		double* phi1, double* phi2, double& scalar){
	__mod_solvers_MOD_product_scalar(&ss1, &ss2, &ss3, &nn1, &nn2, &nn3,
                            phi1,phi2, &scalar);
}
inline void mpt_multadd1(int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3, 
		double& mult, double* vec, double* phi){
	__mod_solvers_MOD_multadd1(&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, &mult, vec, phi);
}
inline void mpt_multadd2(int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3, 
		double& mult1, double& mult2, double* vec, double* phi1, double* phi2, int& init_yes){
	__mod_solvers_MOD_multadd2(&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, 
			&mult1, &mult2, vec, phi1, phi2, &init_yes);
}
inline void mpt_status_iteration(double& eps, double& norm, int& counter, int& n_restarts, int& exit_yes,
		int& quiet_yes1, int& quiet_yes2){
	__mod_solvers_MOD_status_iteration(&eps, &norm, &counter, &n_restarts, &exit_yes, &quiet_yes1, &quiet_yes2);
}
inline void mpt_handle_corner_rhs(int& g, double* phi){__mod_solvers_MOD_handle_corner_rhs(&g, phi);}
inline void mpt_bicgstab2(double& eps, int& n_it_max, int& init_yes, int& n1, int& n2, int& n3, int& g,
		int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3,
		double* bb, double* phi, int& problem_type, int& quiet_yes1, int& quiet_yes2, int& preconditioner){
	__mod_solvers_MOD_bicgstab2(&eps, &n_it_max, &init_yes, &n1, &n2, &n3, &g,
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3, bb, phi, &problem_type, &quiet_yes1,
			&quiet_yes2, &preconditioner);
}
inline void mpt_get_norms2(int& g, int& n1, int& n2, int& n3, 
		int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3,
		double* phi, int& problem_type, int& inf_yes, int& two_yes,
		double& norminf, double& normtwo){
	__mod_solvers_MOD_get_norms2(&g, &n1, &n2, &n3, 
			&ss1, &ss2, &ss3, &nn1, &nn2, &nn3,
			phi, &problem_type, &inf_yes, &two_yes, &norminf, &normtwo);
}
inline void mpt_product_scalar2(int& g, int& ss1, int& ss2, int& ss3, int& nn1, int& nn2, int& nn3, 
		double* phi1, double* phi2, double& scalar){
	__mod_solvers_MOD_product_scalar2(&g, &ss1, &ss2, &ss3, &nn1, &nn2, &nn3,
			phi1, phi2, &scalar);
}
inline void mpt_apply_nullspace2(int& g, double* psi, double* res, int& problem_type){
	__mod_solvers_MOD_apply_nullspace2(&g, psi, res, &problem_type);
}
inline void mpt_relax_restrict(int& init_yes, int& nullspace_yes, int& g,
		double* psi, double* bb, double* phi, double* work, double* coarse, int& problem_type){
	__mod_solvers_MOD_relax_restrict(&init_yes, &nullspace_yes, &g, psi, bb, phi, work, coarse, &problem_type); 
}
inline void mpt_interpolate_relax(int& g, double* bb, double* phi, double* work, double* coarse, int& problem_type){
	__mod_solvers_MOD_interpolate_relax(&g, bb, phi, work, coarse, &problem_type);
}
inline void mpt_plain_restrict(int& nullspace_yes, int& g, 
		double* psi, double* bb, double* work, double* coarse, int& problem_type){
	__mod_solvers_MOD_plain_restrict(&nullspace_yes, &g, psi, bb, work, coarse, &problem_type);
}
inline void mpt_interpolate_mg(int& g, double* bb, double* phi, double* work, double* coarse, int& problem_type){
	__mod_solvers_MOD_interpolate_mg(&g, bb, phi, work, coarse, &problem_type);
}
inline void mpt_relax_bottom(int& init_yes, int& nullspace_yes, int& g, 
		double* psi, double* bb, double* phi, int& problem_type){
	__mod_solvers_MOD_relax_bottom(&init_yes, &nullspace_yes, &g, psi, bb, phi, &problem_type);
}
#endif //MOD_SOLVERS_WRAPPER_H
