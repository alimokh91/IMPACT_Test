#ifndef USR_FUNC_WRAPPER_H
#define USR_FUNC_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	double __usr_func_MOD_smooth_step(int* step_yes, double* xx);
	double __usr_func_MOD_linear_step(double* xx);
	void __usr_func_MOD_fringe_coeff(double* lam_max, double* x_start, double* x_end,
			double* d_rise, double* d_fall, double* xx, double* lam_fringe);
	void __usr_func_MOD_poiseuille_parabola(double* x_start, double* x_end, double* xx,
			double* parab_val);
	double __usr_func_MOD_interface(double* xx);
	double __usr_func_MOD_erf(double* x);
	double __usr_func_MOD_atanh(double* x);
	void __usr_func_MOD_coord_tanh(double* lmax, double* iimax, double* ii0l,
			double* ii0u, double* ii, double* xx, double* dx);
	void __usr_func_MOD_init_hdf5(void);
	void __usr_func_MOD_init_mpi(void);
	void __usr_func_MOD_finl_hdf5(void);
	void __usr_func_MOD_finl_mpi(void);
	void __usr_func_MOD_print_fcomm_size(void);
	void __usr_func_MOD_mpi_bcast_fort(void);
	void __usr_func_MOD_windkessel_integration_step(void);
	void __usr_func_MOD_pdot3EWK(double* pb, double* dpb_dt);
	void __usr_func_MOD_pdot4EWK(double* p, double* dp_dt);
	double __usr_func_MOD_flow_and_rate(int* dir, double* pos, double* center,
			double* radius);
	double __usr_func_MOD_flow_and_rate2D(int* dir, double* pos, double* center,
			double* radius);
	double __usr_func_MOD_flow_and_rate3D(int* dir, double* pos, double* center,
			double* radius);
	double __usr_func_MOD_fitted_flow(double* xx);
	void __usr_func_MOD_apply_fringe_forcing(void);
	void __usr_func_MOD_apply_windkessel_loading(void);
  int __usr_func_MOD_total_n_local_tet_elements(void);
  int __usr_func_MOD_total_n_local_tri_elements(void);
  int __usr_func_MOD_total_n_local_quad_elements(void);
  int __usr_func_MOD_total_n_local_hex_elements(void);
  void __usr_func_MOD_local_to_global_node_id_loc_con_allperiodic2(int* loc_node_id,
      int* glob_node_id);
  void __usr_func_MOD_global_to_local_node_id_loc_con_allperiodic2(int* loc_node_id,
      int* glob_node_id);
  void __usr_func_MOD_local_to_global_node_id_loc_con(int* loc_node_id,
      int* glob_node_id);
  void __usr_func_MOD_global_to_local_node_id_loc_con(int* loc_node_id,
      int* glob_node_id);
  void __usr_func_MOD_global_id2cart_loc_con(int* node_id, int* i, int* j, int* k);
  void __usr_func_MOD_global_cart2id_loc_con(int* node_id, int* i, int* j, int* k);
  int __usr_func_MOD_n_local_nodes(void);
  int __usr_func_MOD_n_global_nodes(void);
  void __usr_func_MOD_local_id2cart_loc_con(int* node_id, int* i, int* j, int* k,
      int* n_local, int* n_ghost);
  void __usr_func_MOD_local_cart2id_loc_con(int* node_id, int* i, int* j, int* k,
      int* n_local, int* n_ghost);
  void __usr_func_MOD_local_tet_element_nodes(int* elem_no, int* node_no, int* node_id);
  void __usr_func_MOD_local_tri_element_nodes(int* elem_no, int* node_no, int* node_id);
  void __usr_func_MOD_local_quad_element_nodes(int* elem_no, int* node_no, int* node_id);
  void __usr_func_MOD_local_hex_element_nodes(int* elem_no, int* node_no, int* node_id);
  void __usr_func_MOD_interpolate_force_pre_vel(int* dir);
  int __usr_func_MOD_block_id(int* ib1, int* ib2, int* ib3);
  void __usr_func_MOD_block_cart(int* block_id, int* ib, int* jb, int* kb);
  int __usr_func_MOD_total_n_global_tet_elements(void);
  int __usr_func_MOD_total_n_global_tri_elements(void);
  int __usr_func_MOD_total_n_global_quad_elements(void);
  int __usr_func_MOD_total_n_global_hex_elements(void);
  int __usr_func_MOD_local_to_global_tet_elem_id(int* loc_elem_id);
  int __usr_func_MOD_local_to_global_tri_elem_id(int* loc_elem_id);
  int __usr_func_MOD_local_to_global_quad_elem_id(int* loc_elem_id);
  int __usr_func_MOD_local_to_global_hex_elem_id(int* loc_elem_id);
  void __usr_func_MOD_residual2volume_force(int* dir);
  void __usr_func_MOD_set_pointers_to_non_allocatable_arrays(void);
  double __usr_func_MOD_rhs_l2_norm(void);
  double __usr_func_MOD_force_l2_norm(void);
  void __usr_func_MOD_compute_and_store_global_block_sizes(int* vel_grid_yes, int* boundary_yes, int* dir);
  void __usr_func_MOD_reset_force_density_component(int* dir);
  void __usr_func_MOD_open_log_iterations_fortran(void);
  void __usr_func_MOD_close_log_iterations_fortran(void);
  void __usr_func_MOD_save_old_velocity(void);
  void __usr_func_MOD_restore_old_velocity(void);

#ifdef __cplusplus
}
#endif

// function wrappers
inline double mpt_smooth_step(int& step_yes, double& xx){
	return __usr_func_MOD_smooth_step(&step_yes, &xx);
}
inline double mpt_linear_step(double& xx){
	return __usr_func_MOD_linear_step(&xx);
}
inline void mpt_fringe_coeff(double& lam_max, double& x_start, double& x_end,
		double& d_rise, double& d_fall, double& xx, double& lam_fringe){
	__usr_func_MOD_fringe_coeff(&lam_max, &x_start, &x_end, &d_rise, &d_fall,
			&xx, &lam_fringe);
}
inline void mpt_poiseuille_parabola(double& x_start, double& x_end, double& xx,
		double& parab_val){
	__usr_func_MOD_poiseuille_parabola(&x_start, &x_end, &xx, &parab_val);
}
inline double mpt_interface(double& xx){return __usr_func_MOD_interface(&xx);}
inline double mpt_erf(double& x){return __usr_func_MOD_erf(&x);}
inline double mpt_atanh(double& x){return __usr_func_MOD_atanh(&x);}
inline void mpt_coord_tanh(double& lmax, double& iimax, double& ii0l, double& ii0u,
		double& ii, double& xx, double& dx){
	__usr_func_MOD_coord_tanh(&lmax, &iimax, &ii0l, &ii0u, &ii, &xx, &dx);
}
inline void mpt_init_hdf5(){__usr_func_MOD_init_hdf5();}
inline void mpt_init_mpi(){__usr_func_MOD_init_mpi();}
inline void mpt_finl_hdf5(){__usr_func_MOD_finl_hdf5();}
inline void mpt_finl_mpi(){__usr_func_MOD_finl_mpi();}
inline void mpt_print_fcomm_size(){__usr_func_MOD_print_fcomm_size();}
inline void mpt_mpi_bcast_fort(){__usr_func_MOD_mpi_bcast_fort();}
inline void mpt_windkessel_integration_step(){__usr_func_MOD_windkessel_integration_step();}
inline void mpt_pdot3EWK(double& pb, double& dpb_dt){__usr_func_MOD_pdot3EWK(&pb, &dpb_dt);}
inline void mpt_pdot4EWK(double* p, double* dp_dt){__usr_func_MOD_pdot4EWK(p, dp_dt);}
inline double mpt_flow_and_rate(int& dir, double& pos, double* center, double& radius){
	return __usr_func_MOD_flow_and_rate(&dir, &pos, center, &radius);
}
inline double mpt_flow_and_rate2D(int& dir, double& pos, double* center, double& radius){
	return __usr_func_MOD_flow_and_rate2D(&dir, &pos, center, &radius);
}
inline double mpt_flow_and_rate3D(int& dir, double& pos, double* center, double& radius){
	return __usr_func_MOD_flow_and_rate3D(&dir, &pos, center, &radius);
}
inline double mpt_fitted_flow(double& xx){return __usr_func_MOD_fitted_flow(&xx);}
inline void mpt_apply_fringe_forcing(){__usr_func_MOD_apply_fringe_forcing();}
inline void mpt_apply_windkessel_loading(){__usr_func_MOD_apply_windkessel_loading();}
inline int mpt_total_n_local_tet_elements(){return __usr_func_MOD_total_n_local_tet_elements();}
inline int mpt_total_n_local_tri_elements(){return __usr_func_MOD_total_n_local_tri_elements();}
inline int mpt_total_n_local_quad_elements(){return __usr_func_MOD_total_n_local_quad_elements();}
inline int mpt_total_n_local_hex_elements(){return __usr_func_MOD_total_n_local_hex_elements();}
inline void mpt_local_to_global_node_id_loc_con_allperiodic2(int& loc_node_id, 
    int& glob_node_id){
  __usr_func_MOD_local_to_global_node_id_loc_con_allperiodic2(&loc_node_id, &glob_node_id);
}
inline void mpt_global_to_local_node_id_loc_con_allperiodic2(int& loc_node_id,
    int& glob_node_id){
  __usr_func_MOD_global_to_local_node_id_loc_con_allperiodic2(&loc_node_id, &glob_node_id);
}
inline void mpt_local_to_global_node_id_loc_con(int& loc_node_id, 
    int& glob_node_id){
  __usr_func_MOD_local_to_global_node_id_loc_con(&loc_node_id, &glob_node_id);
}
inline void mpt_global_to_local_node_id_loc_con(int& loc_node_id,
    int& glob_node_id){
  __usr_func_MOD_global_to_local_node_id_loc_con(&loc_node_id, &glob_node_id);
}
inline void mpt_global_id2cart_loc_con(int& node_id, int& i, int& j, int& k){
  __usr_func_MOD_global_id2cart_loc_con(&node_id, &i, &j, &k);
}
inline void mpt_global_cart2id_loc_con(int& node_id, int& i, int& j, int& k){
  __usr_func_MOD_global_cart2id_loc_con(&node_id, &i, &j, &k);
}
inline int mpt_n_local_nodes(){return __usr_func_MOD_n_local_nodes();}
inline int mpt_n_global_nodes(){return __usr_func_MOD_n_global_nodes();}
inline void mpt_local_id2cart_loc_con(int& node_id, int& i, int& j, int& k,
    int& n_local, int& n_ghost){
  __usr_func_MOD_local_id2cart_loc_con(&node_id, &i, &j, &k, &n_local, &n_ghost);
}
inline void mpt_local_cart2id_loc_con(int& node_id, int& i, int& j, int& k,
    int& n_local, int& n_ghost){
  __usr_func_MOD_local_cart2id_loc_con(&node_id, &i, &j, &k, &n_local, &n_ghost);
}
inline void mpt_local_tet_element_nodes(int& elem_no, int& node_no, int& node_id){
  __usr_func_MOD_local_tet_element_nodes(&elem_no, &node_no, &node_id);
}
inline void mpt_local_tri_element_nodes(int& elem_no, int& node_no, int& node_id){
  __usr_func_MOD_local_tri_element_nodes(&elem_no, &node_no, &node_id);
}
inline void mpt_local_quad_element_nodes(int& elem_no, int& node_no, int& node_id){
  __usr_func_MOD_local_quad_element_nodes(&elem_no, &node_no, &node_id);
}
inline void mpt_local_hex_element_nodes(int& elem_no, int& node_no, int& node_id){
  __usr_func_MOD_local_hex_element_nodes(&elem_no, &node_no, &node_id);
}
inline void mpt_interpolate_force_pre_vel(int& dir){__usr_func_MOD_interpolate_force_pre_vel(&dir);}
inline int mpt_block_id(int& ib1, int& ib2, int& ib3){return __usr_func_MOD_block_id(&ib1, &ib2, &ib3);}
inline void mpt_block_cart(int& block_id, int& ib, int& jb, int& kb){
  __usr_func_MOD_block_cart(&block_id, &ib, &jb, &kb);
}
inline int mpt_total_n_global_tet_elements(){return __usr_func_MOD_total_n_global_tet_elements();}
inline int mpt_total_n_global_tri_elements(){return __usr_func_MOD_total_n_global_tri_elements();}
inline int mpt_total_n_global_quad_elements(){return __usr_func_MOD_total_n_global_quad_elements();}
inline int mpt_total_n_global_hex_elements(){return __usr_func_MOD_total_n_global_hex_elements();}
inline int mpt_local_to_global_tet_elem_id(int& loc_elem_id){
  return __usr_func_MOD_local_to_global_tet_elem_id(&loc_elem_id);
}
inline int mpt_local_to_global_tri_elem_id(int& loc_elem_id){
  return __usr_func_MOD_local_to_global_tri_elem_id(&loc_elem_id);
}
inline int mpt_local_to_global_quad_elem_id(int& loc_elem_id){
  return __usr_func_MOD_local_to_global_quad_elem_id(&loc_elem_id);
}
inline int mpt_local_to_global_hex_elem_id(int& loc_elem_id){
  return __usr_func_MOD_local_to_global_hex_elem_id(&loc_elem_id);
}
inline void mpt_residual2volume_force(int& dir){__usr_func_MOD_residual2volume_force(&dir);}
inline void mpt_set_pointers_to_non_allocatable_arrays(){
  __usr_func_MOD_set_pointers_to_non_allocatable_arrays();
}
inline double mpt_rhs_l2_norm(){return __usr_func_MOD_rhs_l2_norm();}
inline double mpt_force_l2_norm(){return __usr_func_MOD_force_l2_norm();}
inline void mpt_compute_and_store_global_block_sizes(int& vel_grid_yes, int& boundary_yes, int& dir){__usr_func_MOD_compute_and_store_global_block_sizes(&vel_grid_yes, &boundary_yes, &dir);}
inline void mpt_reset_force_density_component(int& dir){
  __usr_func_MOD_reset_force_density_component(&dir);
}
inline void mpt_open_log_iterations_fortran(){
  __usr_func_MOD_open_log_iterations_fortran();
}
inline void mpt_close_log_iterations_fortran(){
  __usr_func_MOD_close_log_iterations_fortran();
}
inline void mpt_save_old_velocity(){
  __usr_func_MOD_save_old_velocity();
}
inline void mpt_restore_old_velocity(){
  __usr_func_MOD_restore_old_velocity();
}

#endif //USR_FUNC_WRAPPER_H
