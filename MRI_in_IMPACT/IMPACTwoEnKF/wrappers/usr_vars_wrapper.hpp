#ifndef USR_VARS_WRAPPER_H
#define USR_VARS_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	extern double __usr_vars_MOD_pi;
	extern int __usr_vars_MOD_wallnorm_dir;
	extern int __usr_vars_MOD_amax;
	extern int __usr_vars_MOD_bmax;
	extern double __usr_vars_MOD_energy_visc;
	extern double __usr_vars_MOD_diss_viscint_old;
	extern double __usr_vars_MOD_y1_origin;
	extern double __usr_vars_MOD_y2_origin;
	extern double __usr_vars_MOD_y3_origin;
        extern int __usr_vars_MOD_elem_type;
	extern int __usr_vars_MOD_fringe_yes;
	extern double __usr_vars_MOD_fringe_start;
	extern double __usr_vars_MOD_fringe_end;
	extern double __usr_vars_MOD_fringe_rise;
	extern double __usr_vars_MOD_fringe_fall;
	extern int __usr_vars_MOD_fringe_dir;
	extern double __usr_vars_MOD_fringe_amp;
	extern double __usr_vars_MOD_fringe_radius;
	extern int __usr_vars_MOD_step_yes;
	extern double __usr_vars_MOD_t_rampup;
        extern int __usr_vars_MOD_parabola_yes;
        extern int __usr_vars_MOD_d3_channel_parab_yes;
        extern int __usr_vars_MOD_pressure_yes;
        extern double __usr_vars_MOD_p_freq;
        extern double __usr_vars_MOD_p_amp;
	extern double __usr_vars_MOD_q_av;
	extern double __usr_vars_MOD_dq_av_dt;
	extern double __usr_vars_MOD_q_av_0;
	extern double __usr_vars_MOD_q_av_1;
	extern double __usr_vars_MOD_q_av_2;
	extern double __usr_vars_MOD_q_av_3;
	extern int __usr_vars_MOD_wk_yes;
	extern int __usr_vars_MOD_wk_type;
	extern int __usr_vars_MOD_wk_flow_dir;
	extern double __usr_vars_MOD_wk_flow_pos;
	extern double __usr_vars_MOD_wk_flow_radius;
	extern double __usr_vars_MOD_r_c;
	extern double __usr_vars_MOD_r_p;
	extern double __usr_vars_MOD_c_art;
	extern double __usr_vars_MOD_l_art;
	extern double __usr_vars_MOD_wk_frge_start;
	extern double __usr_vars_MOD_wk_frge_end;
	extern double __usr_vars_MOD_wk_frge_rise;
	extern double __usr_vars_MOD_wk_frge_fall;
	extern double __usr_vars_MOD_wk_frge_amp;
        extern double* _fringe_center_c;
        extern double* _wk_flow_center_c;
        extern double* _wk_pre_c;
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
#define mpt_pi __usr_vars_MOD_pi
#define mpt_wallnorm_dir __usr_vars_MOD_wallnorm_dir
#define mpt_amax __usr_vars_MOD_amax
#define mpt_bmax __usr_vars_MOD_bmax
#define mpt_energy_visc __usr_vars_MOD_energy_visc
#define mpt_diss_viscint_old __usr_vars_MOD_diss_viscint_old
#define mpt_y1_origin __usr_vars_MOD_y1_origin
#define mpt_y2_origin __usr_vars_MOD_y2_origin
#define mpt_y3_origin __usr_vars_MOD_y3_origin
#define mpt_elem_type __usr_vars_MOD_elem_type
#define mpt_fringe_yes __usr_vars_MOD_fringe_yes
#define mpt_fringe_start __usr_vars_MOD_fringe_start
#define mpt_fringe_end __usr_vars_MOD_fringe_end
#define mpt_fringe_rise __usr_vars_MOD_fringe_rise
#define mpt_fringe_fall __usr_vars_MOD_fringe_fall
#define mpt_fringe_dir __usr_vars_MOD_fringe_dir
#define mpt_fringe_amp __usr_vars_MOD_fringe_amp
#define mpt_fringe_radius __usr_vars_MOD_fringe_radius
#define mpt_step_yes __usr_vars_MOD_step_yes
#define mpt_t_rampup __usr_vars_MOD_t_rampup
#define mpt_parabola_yes __usr_vars_MOD_parabola_yes
#define mpt_3d_channel_parab_yes __usr_vars_MOD_d3_channel_parab_yes
#define mpt_pressure_yes __usr_vars_MOD_pressure_yes
#define mpt_p_freq __usr_vars_MOD_p_freq
#define mpt_p_amp __usr_vars_MOD_p_amp
#define mpt_q_av __usr_vars_MOD_q_av
#define mpt_dq_av_dt __usr_vars_MOD_dq_av_dt
#define mpt_q_av_0 __usr_vars_MOD_q_av_0
#define mpt_q_av_1 __usr_vars_MOD_q_av_1
#define mpt_q_av_2 __usr_vars_MOD_q_av_2
#define mpt_q_av_3 __usr_vars_MOD_q_av_3
#define mpt_wk_yes __usr_vars_MOD_wk_yes
#define mpt_wk_type __usr_vars_MOD_wk_type
#define mpt_wk_flow_dir __usr_vars_MOD_wk_flow_dir
#define mpt_wk_flow_pos __usr_vars_MOD_wk_flow_pos
#define mpt_wk_flow_radius __usr_vars_MOD_wk_flow_radius
#define mpt_r_c __usr_vars_MOD_r_c
#define mpt_r_p __usr_vars_MOD_r_p
#define mpt_c_art __usr_vars_MOD_c_art
#define mpt_l_art __usr_vars_MOD_l_art
#define mpt_wk_frge_start __usr_vars_MOD_wk_frge_start
#define mpt_wk_frge_end __usr_vars_MOD_wk_frge_end
#define mpt_wk_frge_rise __usr_vars_MOD_wk_frge_rise
#define mpt_wk_frge_fall __usr_vars_MOD_wk_frge_fall
#define mpt_wk_frge_amp __usr_vars_MOD_wk_frge_amp
#define mpt_fringe_center(i) _fringe_center_c[i - 1]
#define mpt_wk_flow_center(i) _wk_flow_center_c[i - 1]
#define mpt_wk_pre(i) _wk_pre_c[i - 1]

#endif //USR_VARS_WRAPPER_H
