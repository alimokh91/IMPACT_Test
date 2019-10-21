#ifndef MOD_VARS_WRAPPER_H
#define MOD_VARS_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	// fortran module variables
	// spatial dimensions
	extern int __mod_vars_MOD_dimens;
	// domain and block specifications
	extern int __mod_vars_MOD_n1;
	extern int __mod_vars_MOD_n2;
	extern int __mod_vars_MOD_n3;
        // number of coarse grids
	extern int _n_grids_max_c;
	extern int __mod_vars_MOD_n_grids;
	extern int __mod_vars_MOD_n_grids_limit;
	// dimensions
	extern int _dim_ncb1c_c;
	extern int _dim_ncb1g_c;
	extern int _dim_ncb1d_c;

	extern int _dim_ncb2c_c;
	extern int _dim_ncb2g_c;
	extern int _dim_ncb2d_c;

	extern int _dim_ncb3c_c;
	extern int _dim_ncb3g_c;
	extern int _dim_ncb3d_c;
	// number of stencil coefficients (field)
	// number of coefficients in the field (central differences assumed)
	extern int _nc1c_c;
	extern int _nc1s_c;
	
        extern int _nc2c_c;
	extern int _nc2s_c;

	extern int _nc3c_c;
	extern int _nc3s_c;
	// interval limits of differential coefficient arrays
	extern int _b1u_c;
	extern int _b2u_c;
	extern int _b3u_c;

	extern int _b1l_c;
	extern int _b2l_c;
	extern int _b3l_c;
	// upwind (non-linear)
	extern int _n1l_c;
	extern int _n2l_c;
	extern int _n3l_c;
	
	extern int _n1u_c;
	extern int _n2u_c;
	extern int _n3u_c;
	// divergence
	extern int _d1l_c;
	extern int _d2l_c;
	extern int _d3l_c;
	
	extern int _d1u_c;
	extern int _d2u_c;
	extern int _d3u_c;
	// gradient
	extern int _g1l_c;
	extern int _g2l_c;
	extern int _g3l_c;
	
	extern int _g1u_c;
	extern int _g2u_c;
	extern int _g3u_c;
	
  // grid specifications
	// physical coordinates (global)
	extern double* _y1p_c;
	extern double* _y1u_c;
	
	extern double* _y2p_c;
	extern double* _y2v_c;
	
	extern double* _y3p_c;
	extern double* _y3w_c;
	// physical coordinates (block)
	extern double* _x1p_c;
	extern double* _x1u_c;

	extern double* _x2p_c;
	extern double* _x2v_c;
	
	extern double* _x3p_c;
	extern double* _x3w_c;
	// grid widths (global)	
	extern double* _dy1p_c;
	extern double* _dy1u_c;
	
	extern double* _dy2p_c;
	extern double* _dy2v_c;
	
	extern double* _dy3p_c;
	extern double* _dy3w_c;
	// grid widths (block)
	extern double* _dx1p_c;
	extern double* _dx1u_c;

        extern double* _dx2p_c;
	extern double* _dx2v_c;

	extern double* _dx3p_c;
	extern double* _dx3w_c;

	// work fields
	// velocities
	extern double* _vel_c;
	// non-linear term
	extern double* _nl_c;
	// right hand side
	extern double* _rhs_c;
	// pressure
	extern double* _pre_c;

	// helper fields (pressure iteration)
	extern double* _work1_c;
	extern double* _work2_c;
	extern double* _work3_c;

        extern int* _outlet_c;
        extern int* _ib_c;

	// indexing (interval limits, shifts)
	// index shift (block -> global)
	extern int __mod_vars_MOD_ishift;
	extern int __mod_vars_MOD_jshift;
	extern int __mod_vars_MOD_kshift;
	// domain size (cleaned from periodicity)
	extern int __mod_vars_MOD_dim1;
	extern int __mod_vars_MOD_dim2;
	extern int __mod_vars_MOD_dim3;
	// pressure / concentrations (incl. border)
	extern int __mod_vars_MOD_s1p;
	extern int __mod_vars_MOD_s2p;
	extern int __mod_vars_MOD_s3p;

	extern int __mod_vars_MOD_n1p;
	extern int __mod_vars_MOD_n2p;
	extern int __mod_vars_MOD_n3p;
	// velocities (incl. border)
	extern int __mod_vars_MOD_s11b;
	extern int __mod_vars_MOD_s21b;
	extern int __mod_vars_MOD_s31b;
	
	extern int __mod_vars_MOD_s12b;
	extern int __mod_vars_MOD_s22b;
	extern int __mod_vars_MOD_s32b;
	
	extern int __mod_vars_MOD_s13b;
	extern int __mod_vars_MOD_s23b;
	extern int __mod_vars_MOD_s33b;

	extern int __mod_vars_MOD_n11b;
	extern int __mod_vars_MOD_n21b;
	extern int __mod_vars_MOD_n31b;
	
	extern int __mod_vars_MOD_n12b;
	extern int __mod_vars_MOD_n22b;
	extern int __mod_vars_MOD_n32b;
	
	extern int __mod_vars_MOD_n13b;
	extern int __mod_vars_MOD_n23b;
	extern int __mod_vars_MOD_n33b;
	// velocities (excl. border)
	extern int __mod_vars_MOD_s11;
	extern int __mod_vars_MOD_s21;
	extern int __mod_vars_MOD_s31;
	
	extern int __mod_vars_MOD_s12;
	extern int __mod_vars_MOD_s22;
	extern int __mod_vars_MOD_s32;
	
	extern int __mod_vars_MOD_s13;
	extern int __mod_vars_MOD_s23;
	extern int __mod_vars_MOD_s33;

	extern int __mod_vars_MOD_n11;
	extern int __mod_vars_MOD_n21;
	extern int __mod_vars_MOD_n31;
	
	extern int __mod_vars_MOD_n12;
	extern int __mod_vars_MOD_n22;
	extern int __mod_vars_MOD_n32;
	
	extern int __mod_vars_MOD_n13;
	extern int __mod_vars_MOD_n23;
	extern int __mod_vars_MOD_n33;
	// coarse grids (multigrid, incl. border)
	extern int __mod_vars_MOD_s1r;
	extern int __mod_vars_MOD_s2r;
	extern int __mod_vars_MOD_s3r;

	extern int __mod_vars_MOD_d1r;
	extern int __mod_vars_MOD_d2r;
	extern int __mod_vars_MOD_d3r;
	// coarse grids (multigrid, excl. border)
	extern int __mod_vars_MOD_s11r;
	extern int __mod_vars_MOD_s22r;
	extern int __mod_vars_MOD_s33r;

	extern int __mod_vars_MOD_d11r;
	extern int __mod_vars_MOD_d22r;
	extern int __mod_vars_MOD_d33r;
	// overlapping convention of the blocks (multigrid, see mod_setup)
	extern int _ls1_c;
	extern int _ls2_c;
	extern int _ls3_c;
	// exchange direction (multigrid)
	extern int __mod_vars_MOD_ex1;
	extern int __mod_vars_MOD_ex2;
	extern int __mod_vars_MOD_ex3;

	// boundary conditions
	// global
	extern int __mod_vars_MOD_bc_1l_global;
	extern int __mod_vars_MOD_bc_1u_global;
	extern int __mod_vars_MOD_bc_2l_global;
	extern int __mod_vars_MOD_bc_2u_global;
	extern int __mod_vars_MOD_bc_3l_global;
	extern int __mod_vars_MOD_bc_3u_global;
	// local (block)
	extern int __mod_vars_MOD_bc_1l;
	extern int __mod_vars_MOD_bc_1u;
	extern int __mod_vars_MOD_bc_2l;
	extern int __mod_vars_MOD_bc_2u;
	extern int __mod_vars_MOD_bc_3l;
	extern int __mod_vars_MOD_bc_3u;
	// physical parameters
	extern double __mod_vars_MOD_l1;
	extern double __mod_vars_MOD_l2;
	extern double __mod_vars_MOD_l3;
	extern double __mod_vars_MOD_l1_half;
	extern double __mod_vars_MOD_l2_half;
	extern double __mod_vars_MOD_l3_half;
	extern double __mod_vars_MOD_l1_amp;
	extern double __mod_vars_MOD_l2_amp;
	extern double __mod_vars_MOD_l3_amp;
	extern double __mod_vars_MOD_re;
	extern double __mod_vars_MOD_acos_yes;
	// numerical parameters
	// general
	extern double __mod_vars_MOD_cfl;
	extern double __mod_vars_MOD_time;
	extern double __mod_vars_MOD_dtime;
	extern double __mod_vars_MOD_subtime;
	extern double __mod_vars_MOD_time_start;
	extern double __mod_vars_MOD_time_end;
	extern double __mod_vars_MOD_dtime_max;
	extern double __mod_vars_MOD_dtime0;
	extern double __mod_vars_MOD_dtime_old;
	extern int __mod_vars_MOD_timestep;
	extern int __mod_vars_MOD_timestep_old;
	extern int __mod_vars_MOD_substep;
	extern int __mod_vars_MOD_n_timesteps;
	extern int __mod_vars_MOD_mapping_yes;
	extern int __mod_vars_MOD_upwind_yes;
	extern int __mod_vars_MOD_euler_yes;
	extern int __mod_vars_MOD_stokes_yes;
	extern int __mod_vars_MOD_twostep_yes;
	extern int _filter_bc_yes_c;
	extern int __mod_vars_MOD_timeint_mode;
	extern int __mod_vars_MOD_forcing_mode;
	extern int __mod_vars_MOD_bulkflow_dir;
	extern int __mod_vars_MOD_n_lp_vel;
	extern int __mod_vars_MOD_n_hp_vel;
	extern double __mod_vars_MOD_chi_vel;
	extern double __mod_vars_MOD_vel_bulk;
	// Runge-Kutta coefficients
        extern double* _ark_c;
        extern double* _brk_c;
	extern int _rk_steps_c;
	// look-up table for stability region of time-integration
	extern int _n_stab_c;
        extern int* _stride_large_c;
        extern int* _stride_med_c;
        extern int* _stride_small_c;
	// Helmholtz pre-factors
	extern double __mod_vars_MOD_thetal;
	extern double __mod_vars_MOD_multl;
	// temporal control
	extern int __mod_vars_MOD_int_dtime;
	extern int __mod_vars_MOD_int_lev_pre;
	
	extern int __mod_vars_MOD_write_large;
	extern int __mod_vars_MOD_write_med;
	extern int __mod_vars_MOD_write_small;
	extern double __mod_vars_MOD_time_out_scal;
	extern double __mod_vars_MOD_dtime_out_scal;
	extern double __mod_vars_MOD_time_out_vect;
	extern double __mod_vars_MOD_dtime_out_vect;
	
	extern int __mod_vars_MOD_write_out_scal;
	extern int __mod_vars_MOD_write_out_vect;
	extern int __mod_vars_MOD_new_dtime;
	extern int __mod_vars_MOD_finish_yes;

	extern int __mod_vars_MOD_write_count;
	extern int __mod_vars_MOD_restart;
        extern char __mod_vars_MOD_restart_char[3];
	extern int __mod_vars_MOD_n_conc_old;
	// further control options
	extern int __mod_vars_MOD_task;
	extern int __mod_vars_MOD_read_nullspace_yes;
	extern int __mod_vars_MOD_nullspace_yes;
	extern int __mod_vars_MOD_nullspace_coarse_yes;
	extern int __mod_vars_MOD_nullspace_ortho_yes;

	extern int __mod_vars_MOD_write_stout_yes;
	extern int __mod_vars_MOD_log_iteration_yes;
	extern int __mod_vars_MOD_write_restart_yes;
	extern int __mod_vars_MOD_write_lambda2_yes;
        extern int __mod_vars_MOD_write_force_yes;
        extern int __mod_vars_MOD_write_xdmf_yes;
        extern int __mod_vars_MOD_scale_output_yes;
	extern int __mod_vars_MOD_write_test_yes;
	// global running indices
	extern int __mod_vars_MOD_direction;
	// explicit treatment of corners with Dirichlet boundary conditions
	extern int _corner_yes_c;
	// system time
	extern int __mod_vars_MOD_elatime;
	extern int __mod_vars_MOD_day;
	extern int __mod_vars_MOD_hour;
	extern int __mod_vars_MOD_minu;
	extern int __mod_vars_MOD_sec;
	extern int __mod_vars_MOD_msec;
	// iteration parameters
	// termination criterion / absolute accuracy of velocities
	extern double __mod_vars_MOD_epsu;
	extern double __mod_vars_MOD_epsu0;
	// smoother
	extern int __mod_vars_MOD_jacobi_yes;
	// number of max. iterations
	extern int __mod_vars_MOD_n_it_outer;
	extern int __mod_vars_MOD_n_it_poisson;
	extern int __mod_vars_MOD_n_it_helmh_vel;
	// preconditioning (multigrid)
	extern int __mod_vars_MOD_precond_outer;
	extern int __mod_vars_MOD_precond_poisson;
	extern int __mod_vars_MOD_precond_helmh_vel;
	// number of smoothings per grid level (multigrid)
	extern int __mod_vars_MOD_n_relax_down;
	extern int __mod_vars_MOD_n_relax_up;
	extern int __mod_vars_MOD_n_relax_bottom;
	
	extern int __mod_vars_MOD_weighting_yes;

	extern int* _init_pre_c;
	extern int* _init_vel_c;
	extern double* _precoffset0_c;
	extern double* _precratio0_c;
	extern int* _impl_dir_c;
	
	// iteration stats
	extern double __mod_vars_MOD_dtime_average;
	extern int __mod_vars_MOD_number_poisson;
	
	// MPI
	// communicators
	extern int __mod_vars_MOD_comm_cart;
	
	extern int __mod_vars_MOD_comm_slice1;
	extern int __mod_vars_MOD_comm_bar1;
	extern int __mod_vars_MOD_comm_slice2;
	extern int __mod_vars_MOD_comm_bar2;
	extern int __mod_vars_MOD_comm_slice3;
	extern int __mod_vars_MOD_comm_bar3;
	// ranks of processes
	extern int __mod_vars_MOD_rank;
	extern int __mod_vars_MOD_rank_bar1;
	extern int __mod_vars_MOD_rank_slice1;
	extern int __mod_vars_MOD_rank_bar2;
	extern int __mod_vars_MOD_rank_slice2;
	extern int __mod_vars_MOD_rank_bar3;
	extern int __mod_vars_MOD_rank_slice3;
	// ranks of neighboring processes (in cartesian grid)
	extern int __mod_vars_MOD_rank1l;
	extern int __mod_vars_MOD_rank1u;
	extern int __mod_vars_MOD_rank2l;
	extern int __mod_vars_MOD_rank2u;
	extern int __mod_vars_MOD_rank3l;
	extern int __mod_vars_MOD_rank3u;
	// error handle
	extern int __mod_vars_MOD_merror;
	// request handles
	extern int __mod_vars_MOD_req1l;
	extern int __mod_vars_MOD_req1u;
	extern int __mod_vars_MOD_req2l;
	extern int __mod_vars_MOD_req2u;
	extern int __mod_vars_MOD_req3l;
	extern int __mod_vars_MOD_req3u;
	// HDF5
	extern int __mod_vars_MOD_herror;
	//*** bbecsek 2015: IMMERSED BOUNDARY ADDITIONS ***
	extern double* _fd_c;
        extern double* _vel_old_c;
	extern double __mod_vars_MOD_mu_fluid;
	extern double __mod_vars_MOD_rho_fluid;
        extern double __mod_vars_MOD_l_ref;
        extern double __mod_vars_MOD_u_ref;
  
  


#ifdef __cplusplus
}
#endif


// variable macros
/* scalar variables are kept lowercase to avoid confusion.
array variables are represented in Fortran syntax, i.e. ARRAY(i,j), but are made sure to be accessed
the way they were stored in Fortran (column-major / first index fastest) and are kept lowercase, too.
The running indices are kept in FORTRAN manner, i.e. they DO NOT run from 0 to [N-1] but may even have 
negative indices!
-> 1D array R^M      : mpt_array(i)       -> __modulename_MOD_array[(i - foffset_i)]
-> 2D array R^MxN    : mpt_array(i,j)     -> __modulename_MOD_array[(i - foffset_i) + (j - foffset_j)*M]
-> 3D array R^MxNxO  : mpt_array(i,j,k)   -> __modulename_MOD_array[(i - foffset_i) + (j - foffset_j)*M + (k - foffset_k)*M*N]
-> 4D array R^MxNxOxP: mpt_array(i,j,k,l) -> __modulename_MOD_array[(i - foffset_i) + (j - foffset_j)*M + (k - foffset_k)*M*N + (l - foffset_l)*M*N*O]
where: foffset_x = fortran_start_index_of_x
All variables are added a leading "mpt_" for identifiyng them as IMPACT variables.
*/

// spatial dimensions
#define mpt_dimens __mod_vars_MOD_dimens
// domain and block specifications
#define mpt_n1 __mod_vars_MOD_n1
#define mpt_n2 __mod_vars_MOD_n2
#define mpt_n3 __mod_vars_MOD_n3
// number of coarse grids
#define mpt_n_grids_max _n_grids_max_c
#define mpt_n_grids __mod_vars_MOD_n_grids
#define mpt_n_grids_limit __mod_vars_MOD_n_grids_limit
// dimensions
#define mpt_dim_ncb1c _dim_ncb1c_c
#define mpt_dim_ncb1g _dim_ncb1g_c
#define mpt_dim_ncb1d _dim_ncb1d_c

#define mpt_dim_ncb2c _dim_ncb2c_c
#define mpt_dim_ncb2g _dim_ncb2g_c
#define mpt_dim_ncb2d _dim_ncb2d_c

#define mpt_dim_ncb3c _dim_ncb3c_c
#define mpt_dim_ncb3g _dim_ncb3g_c
#define mpt_dim_ncb3d _dim_ncb3d_c
// number of stencil coefficients (field)
// number of coefficients in the field (central differences assumed)
#define mpt_nc1c _nc1c_c
#define mpt_nc1s _nc1s_c
	
#define mpt_nc2c _nc2c_c
#define mpt_nc2s _nc2s_c

#define mpt_nc3c _nc3c_c
#define mpt_nc3s _nc3s_c
// interval limits of differential coefficient arrays
#define mpt_b1u _b1u_c
#define mpt_b2u _b2u_c
#define	mpt_b3u _b3u_c

#define mpt_b1l _b1l_c
#define mpt_b2l _b2l_c
#define mpt_b3l _b3l_c
// upwind (non-linear)
#define mpt_n1l _n1l_c
#define mpt_n2l _n2l_c
#define mpt_n3l _n3l_c
	
#define mpt_n1u _n1u_c
#define mpt_n2u _n2u_c
#define mpt_n3u _n3u_c
// divergence
#define mpt_d1l _d1l_c
#define mpt_d2l _d2l_c
#define mpt_d3l _d3l_c
	
#define mpt_d1u _d1u_c
#define mpt_d2u _d2u_c
#define mpt_d3u _d3u_c
// gradient
#define mpt_g1l _g1l_c
#define mpt_g2l _g2l_c
#define mpt_g3l _g3l_c
	
#define mpt_g1u _g1u_c
#define mpt_g2u _g2u_c
#define mpt_g3u _g3u_c

// grid specifications
// physical coordinates (global)
#define mpt_y1p(i) _y1p_c[i - 1]
#define mpt_y1u(i) _y1u_c[i]
	
#define mpt_y2p(i) _y2p_c[i - 1]
#define mpt_y2v(i) _y2v_c[i]
	
#define mpt_y3p(i) _y3p_c[i - 1]
#define mpt_y3w(i) _y3w_c[i]
// physical coordinates (block)
#define mpt_x1p(i) _x1p_c[i - _b1l_c]
#define mpt_x1u(i) _x1u_c[i - _b1l_c]

#define mpt_x2p(i) _x2p_c[i - _b2l_c]
#define mpt_x2v(i) _x2v_c[i - _b2l_c]
	
#define mpt_x3p(i) _x3p_c[i - _b3l_c]
#define mpt_x3w(i) _x3w_c[i - _b3l_c]
// grid widths (global)	
#define mpt_dy1p(i) _dy1p_c[i - 1]
#define mpt_dy1u(i) _dy1u_c[i]
	
#define mpt_dy2p(i) _dy2p_c[i - 1]
#define mpt_dy2v(i) _dy2v_c[i]
	
#define mpt_dy3p(i) _dy3p_c[i - 1]
#define mpt_dy3w(i) _dy3w_c[i]
// grid widths (block)
#define mpt_dx1p(i) _dx1p_c[i - 1]
#define mpt_dx1u(i) _dx1u_c[i]

#define mpt_dx2p(i) _dx2p_c[i - 1]
#define mpt_dx2v(i) _dx2v_c[i]

#define mpt_dx3p(i) _dx3p_c[i - 1]
#define mpt_dx3w(i) _dx3w_c[i]

// work fields
// velocities
#define mpt_vel(i,j,k,l) _vel_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1) + (l - 1)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)*(__mod_vars_MOD_n3+_b3u_c-_b3l_c+1)]
// non-linear term
#define mpt_nl(i,j,k,l) _nl_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1) + (l - 1)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)*(__mod_vars_MOD_n3+_b3u_c-_b3l_c+1)]
// right hand side
#define mpt_rhs(i,j,k,l) _rhs_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1) + (l - 1)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)*(__mod_vars_MOD_n3+_b3u_c-_b3l_c+1)]
// pressure
#define mpt_pre(i,j,k) _pre_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)]

// helper fields (pressure iteration)
#define mpt_work1(i,j,k) _work1_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)]

#define mpt_work2(i,j,k) _work2_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)]

#define mpt_work3(i,j,k) _work3_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)]

#define mpt_outlet(i,j,k) _outlet_c[(i - 1) + (j - 1)*3 + (k - 1)*2*3 ]

#define mpt_ib(i,j) _ib_c[(i - 1) + (j - 1)*3]

// indexing (interval limits, shifts)
// index shift (block -> global)
#define mpt_ishift __mod_vars_MOD_ishift
#define mpt_jshift __mod_vars_MOD_jshift
#define mpt_kshift __mod_vars_MOD_kshift
// domain size (cleaned from periodicity)
#define mpt_dim1 __mod_vars_MOD_dim1
#define mpt_dim2 __mod_vars_MOD_dim2
#define mpt_dim3 __mod_vars_MOD_dim3
// pressure / concentrations (incl. border)
#define mpt_s1p __mod_vars_MOD_s1p
#define mpt_s2p __mod_vars_MOD_s2p
#define mpt_s3p __mod_vars_MOD_s3p

#define mpt_n1p __mod_vars_MOD_n1p
#define mpt_n2p __mod_vars_MOD_n2p
#define mpt_n3p __mod_vars_MOD_n3p
// velocities (incl. border)
#define mpt_s11b __mod_vars_MOD_s11b
#define mpt_s21b __mod_vars_MOD_s21b
#define mpt_s31b __mod_vars_MOD_s31b
	
#define mpt_s12b __mod_vars_MOD_s12b
#define mpt_s22b __mod_vars_MOD_s22b
#define mpt_s32b __mod_vars_MOD_s32b
	
#define mpt_s13b __mod_vars_MOD_s13b
#define mpt_s23b __mod_vars_MOD_s23b
#define mpt_s33b __mod_vars_MOD_s33b

#define mpt_n11b __mod_vars_MOD_n11b
#define mpt_n21b __mod_vars_MOD_n21b
#define mpt_n31b __mod_vars_MOD_n31b
	
#define mpt_n12b __mod_vars_MOD_n12b
#define mpt_n22b __mod_vars_MOD_n22b
#define mpt_n32b __mod_vars_MOD_n32b
	
#define mpt_n13b __mod_vars_MOD_n13b
#define mpt_n23b __mod_vars_MOD_n23b
#define mpt_n33b __mod_vars_MOD_n33b
// velocities (excl. border)
#define mpt_s11 __mod_vars_MOD_s11
#define mpt_s21 __mod_vars_MOD_s21
#define mpt_s31 __mod_vars_MOD_s31
	
#define mpt_s12 __mod_vars_MOD_s12
#define mpt_s22 __mod_vars_MOD_s22
#define mpt_s32 __mod_vars_MOD_s32
	
#define mpt_s13 __mod_vars_MOD_s13
#define mpt_s23 __mod_vars_MOD_s23
#define mpt_s33 __mod_vars_MOD_s33

#define mpt_n11 __mod_vars_MOD_n11
#define mpt_n21 __mod_vars_MOD_n21
#define mpt_n31 __mod_vars_MOD_n31
	
#define mpt_n12 __mod_vars_MOD_n12
#define mpt_n22 __mod_vars_MOD_n22
#define mpt_n32 __mod_vars_MOD_n32
	
#define mpt_n13 __mod_vars_MOD_n13
#define mpt_n23 __mod_vars_MOD_n23
#define mpt_n33 __mod_vars_MOD_n33
// coarse grids (multigrid, incl. border)
#define mpt_s1r __mod_vars_MOD_s1r
#define mpt_s2r __mod_vars_MOD_s2r
#define mpt_s3r __mod_vars_MOD_s3r

#define mpt_d1r __mod_vars_MOD_d1r
#define mpt_d2r __mod_vars_MOD_d2r
#define mpt_d3r __mod_vars_MOD_d3r
// coarse grids (multigrid, excl. border)
#define mpt_s11r __mod_vars_MOD_s11r
#define mpt_s22r __mod_vars_MOD_s22r
#define mpt_s33r __mod_vars_MOD_s33r

#define mpt_d11r __mod_vars_MOD_d11r
#define mpt_d22r __mod_vars_MOD_d22r
#define mpt_d33r __mod_vars_MOD_d33r
// overlapping convention of the blocks (multigrid, see mod_setup)
#define mpt_ls1 _ls1_c
#define mpt_ls2 _ls2_c
#define mpt_ls3 _ls3_c
// exchange direction (multigrid)
#define mpt_ex1 __mod_vars_MOD_ex1
#define mpt_ex2 __mod_vars_MOD_ex2
#define mpt_ex3 __mod_vars_MOD_ex3

// boundary conditions
// global
#define mpt_bc_1l_global __mod_vars_MOD_bc_1l_global
#define mpt_bc_1u_global __mod_vars_MOD_bc_1u_global
#define mpt_bc_2l_global __mod_vars_MOD_bc_2l_global
#define mpt_bc_2u_global __mod_vars_MOD_bc_2u_global
#define mpt_bc_3l_global __mod_vars_MOD_bc_3l_global
#define mpt_bc_3u_global __mod_vars_MOD_bc_3u_global
// local (block)
#define mpt_bc_1l __mod_vars_MOD_bc_1l
#define mpt_bc_1u __mod_vars_MOD_bc_1u
#define mpt_bc_2l __mod_vars_MOD_bc_2l
#define mpt_bc_2u __mod_vars_MOD_bc_2u
#define mpt_bc_3l __mod_vars_MOD_bc_3l
#define mpt_bc_3u __mod_vars_MOD_bc_3u
// physical parameters
#define mpt_l1 __mod_vars_MOD_l1
#define mpt_l2 __mod_vars_MOD_l2
#define mpt_l3 __mod_vars_MOD_l3
#define mpt_l1_half __mod_vars_MOD_l1_half
#define mpt_l2_half __mod_vars_MOD_l2_half
#define mpt_l3_half __mod_vars_MOD_l3_half
#define mpt_l1_amp __mod_vars_MOD_l1_amp
#define mpt_l2_amp __mod_vars_MOD_l2_amp
#define mpt_l3_amp __mod_vars_MOD_l3_amp
#define mpt_re __mod_vars_MOD_re
#define mpt_acos_yes __mod_vars_MOD_acos_yes
// numerical parameters
// general
#define mpt_cfl __mod_vars_MOD_cfl
#define mpt_time __mod_vars_MOD_time
#define mpt_dtime __mod_vars_MOD_dtime
#define mpt_subtime __mod_vars_MOD_subtime
#define mpt_time_start __mod_vars_MOD_time_start
#define mpt_time_end __mod_vars_MOD_time_end
#define mpt_dtime_max __mod_vars_MOD_dtime_max
#define mpt_dtime0 __mod_vars_MOD_dtime0
#define mpt_dtime_old __mod_vars_MOD_dtime_old
#define mpt_timestep __mod_vars_MOD_timestep
#define mpt_timestep_old __mod_vars_MOD_timestep_old
#define mpt_substep __mod_vars_MOD_substep
#define mpt_n_timesteps __mod_vars_MOD_n_timesteps
#define mpt_mapping_yes __mod_vars_MOD_mapping_yes
#define mpt_upwind_yes __mod_vars_MOD_upwind_yes
#define mpt_euler_yes __mod_vars_MOD_euler_yes
#define mpt_stokes_yes __mod_vars_MOD_stokes_yes
#define mpt_twostep_yes __mod_vars_MOD_twostep_yes
#define mpt_filter_bc_yes _filter_bc_yes_c
#define mpt_timeint_mode __mod_vars_MOD_timeint_mode
#define mpt_forcing_mode __mod_vars_MOD_forcing_mode
#define mpt_bulkflow_dir __mod_vars_MOD_bulkflow_dir
#define mpt_n_lp_vel __mod_vars_MOD_n_lp_vel
#define mpt_n_hp_vel __mod_vars_MOD_n_hp_vel
#define mpt_chi_vel __mod_vars_MOD_chi_vel
#define mpt_vel_bulk __mod_vars_MOD_vel_bulk
// Runge-Kutta coefficients
#define mpt_ark(i) _ark_c[i - 1]
#define mpt_brk(i) _brk_c[i - 1]
#define mpt_rk_steps _rk_steps_c
// look-up table for stability region of time-integration
#define mpt_n_stab _n_stab_c
#define mpt_stride_large(i) _stride_large_c[i - 1]
#define mpt_stride_med(i) _stride_med_c[i - 1]
#define mpt_stride_small(i) _stride_small_c[i - 1]
// Helmholtz pre-factors
#define mpt_thetal __mod_vars_MOD_thetal
#define mpt_multl __mod_vars_MOD_multl
// temporal control
#define mpt_int_dtime __mod_vars_MOD_int_dtime
#define mpt_int_lev_pre __mod_vars_MOD_int_lev_pre
	
#define mpt_write_large __mod_vars_MOD_write_large
#define mpt_write_med __mod_vars_MOD_write_med
#define mpt_write_small __mod_vars_MOD_write_small
#define mpt_time_out_scal __mod_vars_MOD_time_out_scal
#define mpt_dtime_out_scal __mod_vars_MOD_dtime_out_scal
#define mpt_time_out_vect __mod_vars_MOD_time_out_vect
#define mpt_dtime_out_vect __mod_vars_MOD_dtime_out_vect
	
#define mpt_write_out_scal __mod_vars_MOD_write_out_scal
#define mpt_write_out_vect __mod_vars_MOD_write_out_vect
#define mpt_new_dtime __mod_vars_MOD_new_dtime
#define mpt_finish_yes __mod_vars_MOD_finish_yes

#define mpt_write_count __mod_vars_MOD_write_count
#define mpt_restart __mod_vars_MOD_restart
#define mpt_restart_char __mod_vars_MOD_restart_char // has length 3
#define mpt_n_conc_old __mod_vars_MOD_n_conc_old
// further control options
#define mpt_task __mod_vars_MOD_task
#define mpt_read_nullspace_yes __mod_vars_MOD_read_nullspace_yes
#define mpt_nullspace_yes __mod_vars_MOD_nullspace_yes
#define mpt_nullspace_coarse_yes __mod_vars_MOD_nullspace_coarse_yes
#define mpt_nullspace_ortho_yes __mod_vars_MOD_nullspace_ortho_yes

#define mpt_write_stout_yes __mod_vars_MOD_write_stout_yes
#define mpt_log_iteration_yes __mod_vars_MOD_log_iteration_yes
#define mpt_write_restart_yes __mod_vars_MOD_write_restart_yes
#define mpt_write_lambda2_yes __mod_vars_MOD_write_lambda2_yes
#define mpt_write_force_yes __mod_vars_MOD_write_force_yes
#define mpt_write_xdmf_yes __mod_vars_MOD_write_xdmf_yes
#define mpt_scale_output_yes __mod_vars_MOD_scale_output_yes
#define mpt_write_test_yes __mod_vars_MOD_write_test_yes
// global running indices
#define mpt_direction __mod_vars_MOD_direction
// explicit treatment of corners with Dirichlet boundary conditions
#define mpt_corner_yes _corner_yes_c
// system time
#define mpt_elatime __mod_vars_MOD_elatime
#define mpt_day __mod_vars_MOD_day
#define mpt_hour __mod_vars_MOD_hour
#define mpt_minu __mod_vars_MOD_minu
#define mpt_sec __mod_vars_MOD_sec
#define mpt_msec __mod_vars_MOD_msec
// iteration parameters
// termination criterion / absolute accuracy of velocities
#define mpt_epsu __mod_vars_MOD_epsu
#define mpt_epsu0 __mod_vars_MOD_epsu0
// smoother
#define mpt_jacobi_yes __mod_vars_MOD_jacobi_yes
// number of max. iterations
#define mpt_n_it_outer __mod_vars_MOD_n_it_outer
#define mpt_n_it_poisson __mod_vars_MOD_n_it_poisson
#define mpt_n_it_helmh_vel __mod_vars_MOD_n_it_helmh_vel
// preconditioning (multigrid)
#define mpt_precond_outer __mod_vars_MOD_precond_outer
#define mpt_precond_poisson __mod_vars_MOD_precond_poisson
#define mpt_precond_helmh_vel __mod_vars_MOD_precond_helmh_vel
// number of smoothings per grid level (multigrid)
#define mpt_n_relax_down __mod_vars_MOD_n_relax_down
#define mpt_n_relax_up __mod_vars_MOD_n_relax_up
#define mpt_n_relax_bottom __mod_vars_MOD_n_relax_bottom
	
#define mpt_weighting_yes __mod_vars_MOD_weighting_yes

// iteration parameters
#define mpt_init_pre(i) _init_pre_c[i - 1]
#define mpt_init_vel(i) _init_vel_c[i - 1]
#define mpt_precoffset0(i) _precoffset0_c[i - 1]
#define mpt_precratio0(i) _precratio0_c[i - 1]
#define mpt_impl_dir(i) _impl_dir_c[i - 1]

// iteration stats
#define mpt_dtime_average __mod_vars_MOD_dtime_average
#define mpt_number_poisson __mod_vars_MOD_number_poisson
	
// MPI
// communicators
#define mpt_comm_cart __mod_vars_MOD_comm_cart
	
#define mpt_comm_slice1 __mod_vars_MOD_comm_slice1
#define mpt_comm_bar1 __mod_vars_MOD_comm_bar1
#define mpt_comm_slice2 __mod_vars_MOD_comm_slice2
#define mpt_comm_bar2 __mod_vars_MOD_comm_bar2
#define mpt_comm_slice3 __mod_vars_MOD_comm_slice3
#define mpt_comm_bar3 __mod_vars_MOD_comm_bar3
// ranks of processes
#define mpt_rank __mod_vars_MOD_rank
#define mpt_rank_bar1 __mod_vars_MOD_rank_bar1
#define mpt_rank_slice1 __mod_vars_MOD_rank_slice1
#define mpt_rank_bar2 __mod_vars_MOD_rank_bar2
#define mpt_rank_slice2 __mod_vars_MOD_rank_slice2
#define mpt_rank_bar3 __mod_vars_MOD_rank_bar3
#define mpt_rank_slice3 __mod_vars_MOD_rank_slice3
// ranks of neighboring processes (in cartesian grid)
#define mpt_rank1l __mod_vars_MOD_rank1l
#define mpt_rank1u __mod_vars_MOD_rank1u
#define mpt_rank2l __mod_vars_MOD_rank2l
#define mpt_rank2u __mod_vars_MOD_rank2u
#define mpt_rank3l __mod_vars_MOD_rank3l
#define mpt_rank3u __mod_vars_MOD_rank3u
// error handle
#define mpt_merror __mod_vars_MOD_merror
// request handles
#define mpt_req1l __mod_vars_MOD_req1l
#define mpt_req1u __mod_vars_MOD_req1u
#define mpt_req2l __mod_vars_MOD_req2l
#define mpt_req2u __mod_vars_MOD_req2u
#define mpt_req3l __mod_vars_MOD_req3l
#define mpt_req3u __mod_vars_MOD_req3u
// HDF5
#define mpt_herror __mod_vars_MOD_herror
//*** bbecsek 2015: IMMERSED BOUNDARY ADDITIONS ***
#define mpt_fd(i,j,k,l) _fd_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1) + (l - 1)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)*(__mod_vars_MOD_n3+_b3u_c-_b3l_c+1)]
#define mpt_vel_old(i,j,k,l) _vel_old_c[(i - _b1l_c) + (j - _b2l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1) + (k - _b3l_c)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1) + (l - 1)*(__mod_vars_MOD_n1+_b1u_c-_b1l_c+1)*(__mod_vars_MOD_n2+_b2u_c-_b2l_c+1)*(__mod_vars_MOD_n3+_b3u_c-_b3l_c+1)]
#define mpt_mu_fluid __mod_vars_MOD_mu_fluid
#define mpt_rho_fluid __mod_vars_MOD_rho_fluid
#define mpt_l_ref __mod_vars_MOD_l_ref
#define mpt_u_ref __mod_vars_MOD_u_ref



#endif //MOD_VARS_WRAPPER_H
