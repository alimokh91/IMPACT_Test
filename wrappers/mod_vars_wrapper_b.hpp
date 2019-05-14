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
	extern int __mod_vars_MOD_n_grids_max;
	extern int __mod_vars_MOD_n_grids;
	extern int __mod_vars_MOD_n_grids_limit;
	// dimensions
	extern int __mod_vars_MOD_dim_ncb1c;
	extern int __mod_vars_MOD_dim_ncb1g;
	extern int __mod_vars_MOD_dim_ncb1d;

	extern int __mod_vars_MOD_dim_ncb2c;
	extern int __mod_vars_MOD_dim_ncb2g;
	extern int __mod_vars_MOD_dim_ncb2d;

	extern int __mod_vars_MOD_dim_ncb3c;
	extern int __mod_vars_MOD_dim_ncb3g;
	extern int __mod_vars_MOD_dim_ncb3d;
	// number of stencil coefficients (field)
	// number of coefficients in the field (central differences assumed)
	extern int __mod_vars_MOD_nc1c;
	extern int __mod_vars_MOD_nc1s;
	
	extern int __mod_vars_MOD_nc2c;
	extern int __mod_vars_MOD_nc2s;

	extern int __mod_vars_MOD_nc3c;
	extern int __mod_vars_MOD_nc3s;
	// interval limits of differential coefficient arrays
	extern int __mod_vars_MOD_b1u;
	extern int __mod_vars_MOD_b2u;
	extern int __mod_vars_MOD_b3u;

	extern int __mod_vars_MOD_b1l;
	extern int __mod_vars_MOD_b2l;
	extern int __mod_vars_MOD_b3l;
	// upwind (non-linear)
	extern int __mod_vars_MOD_n1l;
	extern int __mod_vars_MOD_n2l;
	extern int __mod_vars_MOD_n3l;
	
	extern int __mod_vars_MOD_n1u;
	extern int __mod_vars_MOD_n2u;
	extern int __mod_vars_MOD_n3u;
	// divergence
	extern int __mod_vars_MOD_d1l;
	extern int __mod_vars_MOD_d2l;
	extern int __mod_vars_MOD_d3l;
	
	extern int __mod_vars_MOD_d1u;
	extern int __mod_vars_MOD_d2u;
	extern int __mod_vars_MOD_d3u;
	// gradient
	extern int __mod_vars_MOD_g1l;
	extern int __mod_vars_MOD_g2l;
	extern int __mod_vars_MOD_g3l;
	
	extern int __mod_vars_MOD_g1u;
	extern int __mod_vars_MOD_g2u;
	extern int __mod_vars_MOD_g3u;
        // differential coefficient arrays
	// 1st derivative (central)
	extern double __mod_vars_MOD_cp1[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cp2[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cp3[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];

	extern double __mod_vars_MOD_cu1[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cv2[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cw3[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];
	// 1st derivative (upwind)
	extern double __mod_vars_MOD_cnp1d[(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cnp2d[(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cnp3d[(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)*(__mod_vars_MOD_n3+1)];

	extern double __mod_vars_MOD_cnp1u[(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cnp2u[(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cnp3u[(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)*(__mod_vars_MOD_n3+1)];

	extern double __mod_vars_MOD_cnu1d[(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cnv2d[(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cnw3d[(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)*(__mod_vars_MOD_n3+1)];

	extern double __mod_vars_MOD_cnu1u[(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cnv2u[(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cnw3u[(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)*(__mod_vars_MOD_n3+1)];
	// divergence
	extern double __mod_vars_MOD_cdu1[(__mod_vars_MOD_d1u-__mod_vars_MOD_d1l+1)*(_mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cdv2[(__mod_vars_MOD_d2u-__mod_vars_MOD_d2l+1)*(_mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cdw3[(__mod_vars_MOD_d3u-__mod_vars_MOD_d3l+1)*(_mod_vars_MOD_n3+1)];
	// divergence (transposed)
	extern double __mod_vars_MOD_cdu1t[(__mod_vars_MOD_g1u-__mod_vars_MOD_g1l+1)*(_mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cdv2t[(__mod_vars_MOD_g2u-__mod_vars_MOD_g2l+1)*(_mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cdw3t[(__mod_vars_MOD_g3u-__mod_vars_MOD_g3l+1)*(_mod_vars_MOD_n3+1)];
	// gradient
	extern double __mod_vars_MOD_cgp1[(__mod_vars_MOD_g1u-__mod_vars_MOD_g1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cgp2[(__mod_vars_MOD_g2u-__mod_vars_MOD_g2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cgp3[(__mod_vars_MOD_g3u-__mod_vars_MOD_g3l+1)*(__mod_vars_MOD_n3+1)];
	// gradient (transposed)
	extern double __mod_vars_MOD_cgp1t[(__mod_vars_MOD_d1u-__mod_vars_MOD_d1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cgp2t[(__mod_vars_MOD_d2u-__mod_vars_MOD_d2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cgp3t[(__mod_vars_MOD_d3u-__mod_vars_MOD_d3l+1)*(__mod_vars_MOD_n3+1)];
	// 2nd derivative (central)
	extern double __mod_vars_MOD_cp11[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cp22[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cp33[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];
	
	extern double __mod_vars_MOD_cu11[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cv22[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cw33[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];
	// interpolation
	extern double __mod_vars_MOD_cipu[(__mod_vars_MOD_g1u-__mod_vars_MOD_g1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cipv[(__mod_vars_MOD_g2u-__mod_vars_MOD_g2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cipw[(__mod_vars_MOD_g3u-__mod_vars_MOD_g3l+1)*(__mod_vars_MOD_n3+1)];

	extern double __mod_vars_MOD_ciup[(__mod_vars_MOD_d1u-__mod_vars_MOD_d1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_civp[(__mod_vars_MOD_d2u-__mod_vars_MOD_d2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_ciwp[(__mod_vars_MOD_d3u-__mod_vars_MOD_d3l+1)*(__mod_vars_MOD_n3+1)];
	// filter
	extern double __mod_vars_MOD_cfp1[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cfp2[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cfp3[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];

	extern double __mod_vars_MOD_cfu1[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cfv2[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cfw3[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];
	// integrator (not for pressure grid)
	extern double __mod_vars_MOD_cint1[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cint2[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cint3[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)*(__mod_vars_MOD_n3+1)];
	// 2nd derivative (multigrid)
	extern double __mod_vars_MOD_cp11r[3*(__mod_vars_MOD_n1+1)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_cp22r[3*(__mod_vars_MOD_n2+1)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_cp33r[3*(__mod_vars_MOD_n3+1)*__mod_vars_MOD_n_grids_max];

	extern double __mod_vars_MOD_cu11r[3*(__mod_vars_MOD_n1+1)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_cv22r[3*(__mod_vars_MOD_n2+1)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_cw33r[3*(__mod_vars_MOD_n3+1)*__mod_vars_MOD_n_grids_max];

	extern double __mod_vars_MOD_cdg1[3*__mod_vars_MOD_n1*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_cdg2[3*__mod_vars_MOD_n2*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_cdg3[3*__mod_vars_MOD_n3*__mod_vars_MOD_n_grids_max];
	// interpolation (multigrid)
	extern double __mod_vars_MOD_ci1[2*__mod_vars_MOD_n1*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_ci2[2*__mod_vars_MOD_n2*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_ci3[2*__mod_vars_MOD_n3*__mod_vars_MOD_n_grids_max];

	extern double __mod_vars_MOD_cih1[2*(__mod_vars_MOD_n1+1)];
	extern double __mod_vars_MOD_cih2[2*(__mod_vars_MOD_n2+1)];
	extern double __mod_vars_MOD_cih3[2*(__mod_vars_MOD_n3+1)];
	// restriction (multigrid)
	extern double __mod_vars_MOD_cr1[3*__mod_vars_MOD_n1*(__mod_vars_MOD_n_grids_max-1)];
	extern double __mod_vars_MOD_cr2[3*__mod_vars_MOD_n2*(__mod_vars_MOD_n_grids_max-1)];
	extern double __mod_vars_MOD_cr3[3*__mod_vars_MOD_n3*(__mod_vars_MOD_n_grids_max-1)];

	extern double __mod_vars_MOD_crest1[(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n1+1)*(__mod_vars_MOD_n_grids_max-1)];
	extern double __mod_vars_MOD_crest2[(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n2+1)*(__mod_vars_MOD_n_grids_max-1)];
	extern double __mod_vars_MOD_crest3[(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*(__mod_vars_MOD_n3+1)*(__mod_vars_MOD_n_grids_max-1)];

	extern double __mod_vars_MOD_crh1[2*__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_crh2[2*__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_crh3[2*__mod_vars_MOD_n3];
	// grid specifications
	// physical coordinates (global)
	extern double __mod_vars_MOD_y1p[__mod_dims_MOD_m1];
	extern double __mod_vars_MOD_y1u[__mod_dims_MOD_m1+1];
	
	extern double __mod_vars_MOD_y2p[__mod_dims_MOD_m2];
	extern double __mod_vars_MOD_y2v[__mod_dims_MOD_m2+1];
	
	extern double __mod_vars_MOD_y3p[__mod_dims_MOD_m3];
	extern double __mod_vars_MOD_y3w[__mod_dims_MOD_m3+1];
	// physical coordinates (block)
	extern double __mod_vars_MOD_x1p[__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l];
	extern double __mod_vars_MOD_x1u[__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l];

	extern double __mod_vars_MOD_x2p[__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l];
	extern double __mod_vars_MOD_x2v[__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l];
	
	extern double __mod_vars_MOD_x3p[__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l];
	extern double __mod_vars_MOD_x3w[__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l];
	// physical coordinates (block, multigrid)
	extern double __mod_vars_MOD_x1pr[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_x1ur[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*__mod_vars_MOD_n_grids_max];
	
	extern double __mod_vars_MOD_x2pr[(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_x2vr[(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*__mod_vars_MOD_n_grids_max];
	
	extern double __mod_vars_MOD_x3pr[(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*__mod_vars_MOD_n_grids_max];
	extern double __mod_vars_MOD_x3wr[(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*__mod_vars_MOD_n_grids_max];
	// grid widths (global)	
	extern double __mod_vars_MOD_dy1p[__mod_dims_MOD_m1];
	extern double __mod_vars_MOD_dy1u[__mod_dims_MOD_m1+1];
	
	extern double __mod_vars_MOD_dy2p[__mod_dims_MOD_m2];
	extern double __mod_vars_MOD_dy2v[__mod_dims_MOD_m2+1];
	
	extern double __mod_vars_MOD_dy3p[__mod_dims_MOD_m3];
	extern double __mod_vars_MOD_dy3w[__mod_dims_MOD_m3+1];
	// grid widths (block)
	extern double __mod_vars_MOD_dx1p[__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_dx1u[__mod_vars_MOD_n1+1];

	extern double __mod_vars_MOD_dx2p[__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_dx2v[__mod_vars_MOD_n2+1];

	extern double __mod_vars_MOD_dx3p[__mod_vars_MOD_n3];
	extern double __mod_vars_MOD_dx3w[__mod_vars_MOD_n3+1];

	extern double __mod_vars_MOD_dx1dm[__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_dx1pm[__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_ddx1pm[__mod_vars_MOD_n1];

	extern double __mod_vars_MOD_dx2dm[__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_dx2pm[__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_ddx2pm[__mod_vars_MOD_n2];

	extern double __mod_vars_MOD_dx3dm[__mod_vars_MOD_n3];
	extern double __mod_vars_MOD_dx3pm[__mod_vars_MOD_n3];
	extern double __mod_vars_MOD_ddx3pm[__mod_vars_MOD_n3];

	extern double __mod_vars_MOD_dx1gm[__mod_vars_MOD_n1+1];
	extern double __mod_vars_MOD_dx1um[__mod_vars_MOD_n1+1];
	extern double __mod_vars_MOD_ddx1um[__mod_vars_MOD_n1+1];

	extern double __mod_vars_MOD_dx2gm[__mod_vars_MOD_n2+1];
	extern double __mod_vars_MOD_dx2vm[__mod_vars_MOD_n2+1];
	extern double __mod_vars_MOD_ddx2vm[__mod_vars_MOD_n2+1];

	extern double __mod_vars_MOD_dx3gm[__mod_vars_MOD_n3+1];
	extern double __mod_vars_MOD_dx3wm[__mod_vars_MOD_n3+1];
	extern double __mod_vars_MOD_ddx3wm[__mod_vars_MOD_n3+1];
	// work fields
	// velocities
	extern double __mod_vars_MOD_vel[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*3];
	// non-linear term
	extern double __mod_vars_MOD_nl[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*3];
	// right hand side
	extern double __mod_vars_MOD_rhs[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*3];
	// pressure
	extern double __mod_vars_MOD_pre[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// outflow boundary conditions (velocity field)
	extern double __mod_vars_MOD_bc11[__mod_vars_MOD_n2*__mod_vars_mod_n3*2];
	extern double __mod_vars_MOD_nlbc11[__mod_vars_MOD_n2*__mod_vars_mod_n3*2];

	extern double __mod_vars_MOD_bc12[(__mod_vars_MOD_n1+1)*__mod_vars_mod_n3*2];
	extern double __mod_vars_MOD_nlbc12[(__mod_vars_MOD_n1+1)*__mod_vars_mod_n3*2];

	extern double __mod_vars_MOD_bc13[(__mod_vars_MOD_n1+1)*__mod_vars_mod_n2*2];
	extern double __mod_vars_MOD_nlbc13[(__mod_vars_MOD_n1+1)*__mod_vars_mod_n2*2];

	extern double __mod_vars_MOD_bc21[(__mod_vars_MOD_n2+1)*__mod_vars_mod_n3*2];
	extern double __mod_vars_MOD_nlbc21[(__mod_vars_MOD_n2+1)*__mod_vars_mod_n3*2];

	extern double __mod_vars_MOD_bc22[__mod_vars_MOD_n1*__mod_vars_mod_n3*2];
	extern double __mod_vars_MOD_nlbc22[__mod_vars_MOD_n1*__mod_vars_mod_n3*2];

	extern double __mod_vars_MOD_bc23[__mod_vars_MOD_n1*(__mod_vars_mod_n2+1)*2];
	extern double __mod_vars_MOD_nlbc23[__mod_vars_MOD_n1*(__mod_vars_mod_n2+1)*2];

	extern double __mod_vars_MOD_bc31[__mod_vars_MOD_n2*(__mod_vars_mod_n3+1)*2];
	extern double __mod_vars_MOD_nlbc31[__mod_vars_MOD_n2*(__mod_vars_mod_n3+1)*2];

	extern double __mod_vars_MOD_bc32[__mod_vars_MOD_n1*(__mod_vars_mod_n3+1)*2];
	extern double __mod_vars_MOD_nlbc32[__mod_vars_MOD_n1*(__mod_vars_mod_n3+1)*2];

	extern double __mod_vars_MOD_bc33[__mod_vars_MOD_n1*__mod_vars_mod_n2*2];
	extern double __mod_vars_MOD_nlbc33[__mod_vars_MOD_n1*__mod_vars_mod_n2*2];

	extern double __mod_vars_MOD_drift1[(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*2];
	extern double __mod_vars_MOD_drift2[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*2];
	extern double __mod_vars_MOD_drift3[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*2];

	// residual
	extern double __mod_vars_MOD_res[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// pressure gradient
	extern double __mod_vars_MOD_gpre[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// weights for divergence-freeness
	extern double __mod_vars_MOD_weight[__mod_vars_MOD_n1*__mod_vars_MOD_n2*__mod_vars_MOD_n3];
	// null-space vector	
	extern double __mod_vars_MOD_psi[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_vel[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*3];

	extern double __mod_vars_MOD_psi_rel1[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// NOTE: for the variables below nn is needed and has not yet been interfaced (may need to be moved above these lines)
	extern double __mod_vars_MOD_psi_rel2[(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*1]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel3[(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*2]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel4[(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*3]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel5[(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*4]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel6[(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*5]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel7[(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*6]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel8[(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*7]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel9[(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*8]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel10[(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*9]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel11[(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*10]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel12[(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*11]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel13[(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*12]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel14[(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*13]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_psi_rel15[(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*14]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*14]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	
	extern double __mod_vars_MOD_th11[__mod_vars_MOD_n2*__mod_vars_MOD_n3*2];
	extern double __mod_vars_MOD_th12[(__mod_vars_MOD_n1+1)*__mod_vars_MOD_n3*2];
	extern double __mod_vars_MOD_th13[(__mod_vars_MOD_n1+1)*__mod_vars_MOD_n2*2];

	extern double __mod_vars_MOD_th21[(__mod_vars_MOD_n2+1)*__mod_vars_MOD_n3*2];
	extern double __mod_vars_MOD_th22[__mod_vars_MOD_n1*__mod_vars_MOD_n3*2];
	extern double __mod_vars_MOD_th23[__mod_vars_MOD_n1*(__mod_vars_MOD_n2+1)*2];

	extern double __mod_vars_MOD_th31[__mod_vars_MOD_n2*(__mod_vars_MOD_n3+1)*2];
	extern double __mod_vars_MOD_th32[__mod_vars_MOD_n1*(__mod_vars_MOD_n3+1)*2];
	extern double __mod_vars_MOD_th33[__mod_vars_MOD_n1*__mod_vars_MOD_n2*2];

	// multigrid
	extern double __mod_vars_MOD_vec1c[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// NOTE: again nn is needed below (move it up?)
	extern double __mod_vars_MOD_vec2a[(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*1]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec2b[(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*1]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec2c[(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*1]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec3a[(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*2]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec3b[(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*2]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec3c[(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*2]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec4a[(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*3]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec4b[(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*3]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec4c[(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*3]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec5a[(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*4]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec5b[(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*4]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec5c[(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*4]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec6a[(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*5]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec6b[(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*5]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec6c[(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*5]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec7a[(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*6]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec7b[(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*6]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec7c[(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*6]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec8a[(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*7]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec8b[(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*7]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec8c[(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*7]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec9a[(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*8]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec9b[(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*8]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec9c[(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*8]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec10a[(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*9]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec10b[(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*9]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec10c[(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*9]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec11a[(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*10]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec11b[(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*10]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec11c[(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*10]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec12a[(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*11]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec12b[(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*11]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec12c[(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*11]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec13a[(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*12]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec13b[(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*12]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec13c[(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*12]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec14a[(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*13]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec14b[(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*13]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec14c[(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*13]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

	extern double __mod_vars_MOD_vec15a[(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*14]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*14]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_vec15b[(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_nn[1+3*14]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_nn[2+3*14]+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// BiCGstab / Richardson
	extern double __mod_vars_MOD_pp[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_ap[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_rr[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_rh[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_ar[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_z1[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_z2[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// product_div_grad
	extern double __mod_vars_MOD_dig[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// helper fields (pressure iteration)
	extern double __mod_vars_MOD_work1[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_work2[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	extern double __mod_vars_MOD_work3[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
	// line relaxation
	extern double __mod_vars_MOD_vec1[__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_dia1[__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_sor1[__mod_vars_MOD_n1];
	extern double __mod_vars_MOD_band1[2*__mod_vars_MOD_n1];

	extern double __mod_vars_MOD_vec2[__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_dia2[__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_sor2[__mod_vars_MOD_n2];
	extern double __mod_vars_MOD_band2[2*__mod_vars_MOD_n2];

	extern double __mod_vars_MOD_vec3[__mod_vars_MOD_n3];
	extern double __mod_vars_MOD_dia3[__mod_vars_MOD_n3];
	extern double __mod_vars_MOD_sor3[__mod_vars_MOD_n3];
	extern double __mod_vars_MOD_band3[2*__mod_vars_MOD_n3];

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
	extern int __mod_vars_MOD_ls1;
	extern int __mod_vars_MOD_ls2;
	extern int __mod_vars_MOD_ls3;
	// exchange direction (multigrid)
	extern int __mod_vars_MOD_ex1;
	extern int __mod_vars_MOD_ex2;
	extern int __mod_vars_MOD_ex3;

	// boundary conditions
	// global
	extern int __mod_vars_MOD_outlet[3*2*3];

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
	// field properties
	extern int __mod_vars_MOD_n_gather[3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_nn[3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_nb[3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_ib[3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_snf[2*3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_snb[2*3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_bc[2*3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_ngb[2*3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_comm1[__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_comm2[__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_rankc2[__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_participate_yes[__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_recvr[__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_recvi[__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_dispr[__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_dispi[__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_offsr[3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_offsi[3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_sizsr[3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	extern int __mod_vars_MOD_sizsi[3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3*__mod_vars_MOD_n_grids_max];
	// physical parameters
	extern double __mod_vars_MOD_l1;
	extern double __mod_vars_MOD_l2;
	extern double __mod_vars_MOD_l3;
	extern double __mod_vars_MOD_re;
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
	extern int __mod_vars_MOD_filter_bc_yes;
	extern int __mod_vars_MOD_timeint_mode;
	extern int __mod_vars_MOD_forcing_mode;
	extern int __mod_vars_MOD_bulkflow_dir;
	extern int __mod_vars_MOD_n_lp_vel;
	extern int __mod_vars_MOD_n_hp_vel;
	extern double __mod_vars_MOD_chi_vel;
	extern double __mod_vars_MOD_vel_bulk;
	// Runge-Kutta coefficients
	extern double __mod_vars_MOD_ark[3];
	extern double __mod_vars_MOD_brk[3];
	extern int __mod_vars_MOD_rk_steps;
	// look-up table for stability region of time-integration
	extern double __mod_vars_MOD_stabilitylimit[41];
	extern int __mod_vars_MOD_n_stab;
	// Helmholtz pre-factors
	extern double __mod_vars_MOD_thetal;
	extern double __mod_vars_MOD_multl;
	// temporal control
	extern int __mod_vars_MOD_int_dtime;
	extern int __mod_vars_MOD_int_lev_pre;
	
	extern int __mod_vars_MOD_stride_large[3];
	extern int __mod_vars_MOD_stride_med[3];
	extern int __mod_vars_MOD_stride_small[3];
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
	extern int __mod_vars_MOD_write_test_yes;
	// global running indices
	extern int __mod_vars_MOD_direction;
	// explicit treatment of corners with Dirichlet boundary conditions
	extern int __mod_vars_MOD_corner_yes;
	// system time
	extern int __mod_vars_MOD_elatime;
	extern int __mod_vars_MOD_ctime[8];
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
	// expected convergence rate (outer iteration)
	extern double __mod_vars_MOD_precratio0[__mod_vars_MOD_rk_steps];
	extern double __mod_vars_MOD_precoffset0[__mod_vars_MOD_rk_steps];
	extern double __mod_vars_MOD_precratio[__mod_vars_MOD_rk_steps*__mod_vars_MOD_n_it_outer];
	extern double __mod_vars_MOD_precoffset[__mod_vars_MOD_rk_steps*__mod_vars_MOD_n_it_outer];
	// null-initialization (outer iteration)
	extern int __mod_vars_MOD_init_pre[__mod_vars_MOD_rk_steps];
	extern int __mod_vars_MOD_init_vel[__mod_vars_MOD_rk_steps];
	// preconditioning (multigrid)
	extern int __mod_vars_MOD_precond_outer;
	extern int __mod_vars_MOD_precond_poisson;
	extern int __mod_vars_MOD_precond_helmh_vel;
	// number of smoothings per grid level (multigrid)
	extern int __mod_vars_MOD_n_relax_down;
	extern int __mod_vars_MOD_n_relax_up;
	extern int __mod_vars_MOD_n_relax_bottom;
	// implicit directions of line relaxations (multigrid)
	extern int __mod_vars_MOD_impl_dir[3];
	
	extern int __mod_vars_MOD_weighting_yes;
	
	// iteration stats
	extern double __mod_vars_MOD_dtime_average;
	extern double __mod_vars_MOD_max_div_init[2];
	extern int __mod_vars_MOD_number_poisson;
	// counters
	extern int __mod_vars_MOD_count0[__mod_vars_MOD_rk_steps];
	extern int __mod_vars_MOD_countp[__mod_vars_MOD_rk_steps*2];
	extern int __mod_vars_MOD_counth[__mod_vars_MOD_rk_steps*3];
	// convergence rate
	extern int __mod_vars_MOD_ratio0[__mod_vars_MOD_rk_steps];
	extern int __mod_vars_MOD_ratiop[__mod_vars_MOD_rk_steps*2];
	extern int __mod_vars_MOD_ratioh[__mod_vars_MOD_rk_steps*3];
	
	// MPI
	// communicators
	extern int __mod_vars_MOD_comm_cart;
	
	extern int __mod_vars_MOD_comm_slice1;
	extern int __mod_vars_MOD_comm_bar1;
	extern int __mod_vars_MOD_comm_slice2;
	extern int __mod_vars_MOD_comm_bar2;
	extern int __mod_vars_MOD_comm_slice3;
	extern int __mod_vars_MOD_comm_bar3;
 	// dimension and position of the blocks inside the communicators (grid indices)
	extern int __mod_vars_MOD_bar1_size[__mod_dims_MOD_nb1];
	extern int __mod_vars_MOD_bar1_offset[__mod_dims_MOD_nb1];
	extern int __mod_vars_MOD_bar2_size[__mod_dims_MOD_nb2];
	extern int __mod_vars_MOD_bar2_offset[__mod_dims_MOD_nb2];
	extern int __mod_vars_MOD_bar3_size[__mod_dims_MOD_nb3];
	extern int __mod_vars_MOD_bar3_offset[__mod_dims_MOD_nb3];
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
	//*** bbecsek 2015: IMMERSED BOUNDARY PARAMETERS ***
	extern double __mod_vars_MOD_reach;
	extern int __mod_vars_MOD_m_bound;
	extern int __mod_vars_MOD_m_elems;
	extern double __mod_vars_MOD_yb[__mod_vars_MOD_m_bound*__mod_vars_MOD_dimens];
	extern double __mod_vars_MOD_xb[__mod_vars_MOD_m_bound*__mod_vars_MOD_dimens];
	extern double __mod_vars_MOD_db[__mod_vars_MOD_m_bound*__mod_vars_MOD_dimens];
	extern double __mod_vars_MOD_ub[__mod_vars_MOD_m_bound*__mod_vars_MOD_dimens];
	extern double __mod_vars_MOD_fb[__mod_vars_MOD_m_bound*__mod_vars_MOD_dimens];
	extern double __mod_vars_MOD_container[__mod_vars_MOD_m_bound*__mod_vars_MOD_dimens];
	extern double __mod_vars_MOD_fd[(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)
		*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)
		*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*3];
 	extern int __mod_vars_MOD_ddf_type;
	extern int __mod_vars_MOD_ib_on;
	extern double __mod_vars_MOD_mu_fluid;
	extern double __mod_vars_MOD_rho_fluid;
	extern double __mod_vars_MOD_e_mod;
	extern double __mod_vars_MOD_nu_poiss;
	extern double __mod_vars_MOD_rho_solid;
	extern double __mod_vars_MOD_strain[__mod_vars_MOD_m_elems*3];
	extern double __mod_vars_MOD_stress[__mod_vars_MOD_m_elems*3];
	extern double __mod_vars_MOD_k_stiff[(2*__mod_vars_MOD_m_bound)*(2*__mod_vars_MOD_m_bound)];
	extern double __mod_vars_MOD_m_mass[(2*__mod_vars_MOD_m_bound)*(2*__mod_vars_MOD_m_bound)];
	extern double __mod_vars_MOD_c_damp[(2*__mod_vars_MOD_m_bound)*(2*__mod_vars_MOD_m_bound)];
	extern double __mod_vars_MOD_c[9];
	extern int __mod_vars_MOD_elems[__mod_vars_MOD_m_elems*3];
	extern double __mod_vars_MOD_l_ref;
	extern double __mod_vars_MOD_u_ref;
	extern int __mod_vars_MOD_req_bcast[3];
	extern int __mod_vars_MOD_req_red[2];
	extern int __mod_vars_MOD_excl_p[3];
	extern double __mod_vars_MOD_delta_b1;
	extern double __mod_vars_MOD_delta_b2;
	extern double __mod_vars_MOD_delta_b3;
	extern int __mod_vars_MOD_total_ranks;
	extern int __mod_vars_MOD_m_ebcs;
	extern double __mod_vars_MOD_ebcs[__mod_vars_MOD_m_ebcs*(__mod_vars_MOD_dimens+1)];
	extern double __mod_vars_MOD_node_vol[__mod_vars_MOD_m_bound];
	extern int __mod_vars_MOD_comm_local;
	extern int __mod_vars_MOD_comm_inter;
	extern double __mod_vars_MOD_max_ev_freq;


#ifdef __cplusplus
}
#endif


// variable macros
/* scalar variables are kept lowercase to avoid confusion.
array variables are represented in Fortran syntax, i.e. ARRAY(i,j), but are made sure to be accessed
the way they were stored in Fortran (column-major / first index fastest) and are kept lowercase, too.
The running indices are kept in FORTRAN manner, i.e. they DO NOT run from 0 to [N-1] but may even have 
negative indices!
-> 1D array R^M      : ARRAY(i)       -> __modulename_MOD_array[(i - foffset_i)]
-> 2D array R^MxN    : ARRAY(i,j)     -> __modulename_MOD_array[(i - foffset_i) + (j - foffset_j)*M]
-> 3D array R^MxNxO  : ARRAY(i,j,k)   -> __modulename_MOD_array[(i - foffset_i) + (j - foffset_j)*M + (k - foffset_k)*M*N]
-> 4D array R^MxNxOxP: ARRAY(i,j,k,l) -> __modulename_MOD_array[(i - foffset_i) + (j - foffset_j)*M + (k - foffset_k)*M*N + (l - foffset_l)*M*N*O]
where: foffset_x = fortran_start_index_of_x
All variables are added a leading "mpt_" for identifiyng them as IMPACT variables.
*/

// spatial dimensions
define mpt_dimens __mod_vars_MOD_dimens;
// domain and block specifications
define mpt_n1 __mod_vars_MOD_n1;
define mpt_n2 __mod_vars_MOD_n2;
define mpt_n3 __mod_vars_MOD_n3;
// number of coarse grids
define mpt_n_grids_max __mod_vars_MOD_n_grids_max;
define mpt_n_grids __mod_vars_MOD_n_grids;
define mpt_n_grids_limit __mod_vars_MOD_n_grids_limit;
// dimensions
define mpt_dim_ncb1c __mod_vars_MOD_dim_ncb1c;
define mpt_dim_ncb1g __mod_vars_MOD_dim_ncb1g;
define mpt_dim_ncb1d __mod_vars_MOD_dim_ncb1d;

define mpt_dim_ncb2c __mod_vars_MOD_dim_ncb2c;
define mpt_dim_ncb2g __mod_vars_MOD_dim_ncb2g;
define mpt_dim_ncb2d __mod_vars_MOD_dim_ncb2d;

define mpt_dim_ncb3c __mod_vars_MOD_dim_ncb3c;
define mpt_dim_ncb3g __mod_vars_MOD_dim_ncb3g;
define mpt_dim_ncb3d __mod_vars_MOD_dim_ncb3d;
// number of stencil coefficients (field)
// number of coefficients in the field (central differences assumed)
define mpt_nc1c __mod_vars_MOD_nc1c;
define mpt_nc1s __mod_vars_MOD_nc1s;
	
define mpt_nc2c __mod_vars_MOD_nc2c;
define mpt_nc2s __mod_vars_MOD_nc2s;

define mpt_nc3c __mod_vars_MOD_nc3c;
define mpt_nc3s __mod_vars_MOD_nc3s;
// interval limits of differential coefficient arrays
define mpt_b1u __mod_vars_MOD_b1u;
define mpt_b2u __mod_vars_MOD_b2u;
define	mpt_b3u __mod_vars_MOD_b3u;

define mpt_b1l __mod_vars_MOD_b1l;
define mpt_b2l __mod_vars_MOD_b2l;
define mpt_b3l __mod_vars_MOD_b3l;
// upwind (non-linear)
define mpt_n1l __mod_vars_MOD_n1l;
define mpt_n2l __mod_vars_MOD_n2l;
define mpt_n3l __mod_vars_MOD_n3l;
	
define mpt_n1u __mod_vars_MOD_n1u;
define mpt_n2u __mod_vars_MOD_n2u;
define mpt_n3u __mod_vars_MOD_n3u;
// divergence
define mpt_d1l __mod_vars_MOD_d1l;
define mpt_d2l __mod_vars_MOD_d2l;
define mpt_d3l __mod_vars_MOD_d3l;
	
define mpt_d1u __mod_vars_MOD_d1u;
define mpt_d2u __mod_vars_MOD_d2u;
define mpt_d3u __mod_vars_MOD_d3u;
// gradient
define mpt_g1l __mod_vars_MOD_g1l;
define mpt_g2l __mod_vars_MOD_g2l;
define mpt_g3l __mod_vars_MOD_g3l;
	
define mpt_g1u __mod_vars_MOD_g1u;
define mpt_g2u __mod_vars_MOD_g2u;
define mpt_g3u __mod_vars_MOD_g3u;
// differential coefficient arrays
// 1st derivative (central)
define mpt_cp1(i,j) __mod_vars_MOD_cp1[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cp2(i,j) __mod_vars_MOD_cp2[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cp3(i,j) __mod_vars_MOD_cp3[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];

define mpt_cu1(i,j) __mod_vars_MOD_cu1[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cv2(i,j) __mod_vars_MOD_cv2[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cv3(i,j) __mod_vars_MOD_cw3[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];
// 1st derivative (upwind)
define mpt_cnp1d(i,j) __mod_vars_MOD_cnp1d[(i - __mod_vars_MOD_n1l) + j*(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)];
define mpt_cnp2d(i,j) __mod_vars_MOD_cnp2d[(i - __mod_vars_MOD_n2l) + j*(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)];
define mpt_cnp3d(i,j) __mod_vars_MOD_cnp3d[(i - __mod_vars_MOD_n3l) + j*(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)];

define mpt_cnp1u(i,j) __mod_vars_MOD_cnp1u[(i - __mod_vars_MOD_n1l) + j*(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)];
define mpt_cnp2u(i,j) __mod_vars_MOD_cnp2u[(i - __mod_vars_MOD_n2l) + j*(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)];
define mpt_cnp3u(i,j) __mod_vars_MOD_cnp3u[(i - __mod_vars_MOD_n3l) + j*(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)];

define mpt_cnu1d(i,j) __mod_vars_MOD_cnu1d[(i - __mod_vars_MOD_n1l) + j*(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)];
define mpt_cnv2d(i,j) __mod_vars_MOD_cnv2d[(i - __mod_vars_MOD_n2l) + j*(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)];
define mpt_cnw3d(i,j) __mod_vars_MOD_cnw3d[(i - __mod_vars_MOD_n3l) + j*(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)];

define mpt_cnu1u(i,j) __mod_vars_MOD_cnu1u[(i - __mod_vars_MOD_n1l) + j*(__mod_vars_MOD_n1u-__mod_vars_MOD_n1l+1)];
define mpt_cnv2u(i,j) __mod_vars_MOD_cnv2u[(i - __mod_vars_MOD_n2l) + j*(__mod_vars_MOD_n2u-__mod_vars_MOD_n2l+1)];
define mpt_cnw3u(i,j) __mod_vars_MOD_cnw3u[(i - __mod_vars_MOD_n3l) + j*(__mod_vars_MOD_n3u-__mod_vars_MOD_n3l+1)];
// divergence
define mpt_cdu1(i,j) __mod_vars_MOD_cdu1[(i - __mod_vars_MOD_d1l) + j*(__mod_vars_MOD_d1u-__mod_vars_MOD_d1l+1)];
define mpt_cdv2(i,j) __mod_vars_MOD_cdv2[(i - __mod_vars_MOD_d2l) + j*(__mod_vars_MOD_d2u-__mod_vars_MOD_d2l+1)];
define mpt_cdw3(i,j) __mod_vars_MOD_cdw3[(i - __mod_vars_MOD_d3l) + j*(__mod_vars_MOD_d3u-__mod_vars_MOD_d3l+1)];
// divergence (transposed)
define mpt_cdu1t(i,j) __mod_vars_MOD_cdu1t[(i - __mod_vars_MOD_g1l) + j*(__mod_vars_MOD_g1u-__mod_vars_MOD_g1l+1)];
define mpt_cdv2t(i,j) __mod_vars_MOD_cdv2t[(i - __mod_vars_MOD_g2l) + j*(__mod_vars_MOD_g2u-__mod_vars_MOD_g2l+1)];
define mpt_cdw3t(i,j) __mod_vars_MOD_cdw3t[(i - __mod_vars_MOD_g3l) + j*(__mod_vars_MOD_g3u-__mod_vars_MOD_g3l+1)];
// gradient
define mpt_cgp1(i,j) __mod_vars_MOD_cgp1[(i - __mod_vars_MOD_g1l) + j*(__mod_vars_MOD_g1u-__mod_vars_MOD_g1l+1)];
define mpt_cgp2(i,j) __mod_vars_MOD_cgp2[(i - __mod_vars_MOD_g2l) + j*(__mod_vars_MOD_g2u-__mod_vars_MOD_g2l+1)];
define mpt_cgp3(i,j) __mod_vars_MOD_cgp3[(i - __mod_vars_MOD_g3l) + j*(__mod_vars_MOD_g3u-__mod_vars_MOD_g3l+1)];
// gradient (transposed)
define mpt_cgp1t(i,j) __mod_vars_MOD_cgp1t[(i - __mod_vars_MOD_d1l) + j*(__mod_vars_MOD_d1u-__mod_vars_MOD_d1l+1)];
define mpt_cgp2t(i,j) __mod_vars_MOD_cgp2t[(i - __mod_vars_MOD_d2l) + j*(__mod_vars_MOD_d2u-__mod_vars_MOD_d2l+1)];
define mpt_cgp2t(i,j) __mod_vars_MOD_cgp3t[(i - __mod_vars_MOD_d3l) + j*(__mod_vars_MOD_d3u-__mod_vars_MOD_d3l+1)];
// 2nd derivative (central)
define mpt_cp11(i,j) __mod_vars_MOD_cp11[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cp22(i,j) __mod_vars_MOD_cp22[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cp33(i,j) __mod_vars_MOD_cp33[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];
	
define mpt_cu11(i,j) __mod_vars_MOD_cu11[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cv22(i,j) __mod_vars_MOD_cv22[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cw33(i,j) __mod_vars_MOD_cw33[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];
// interpolation
define mpt_cipu(i,j) __mod_vars_MOD_cipu[(i - __mod_vars_MOD_g1l) + j*(__mod_vars_MOD_g1u-__mod_vars_MOD_g1l+1)];
define mpt_cipv(i,j) __mod_vars_MOD_cipv[(i - __mod_vars_MOD_g2l) + j*(__mod_vars_MOD_g2u-__mod_vars_MOD_g2l+1)];
define mpt_cipw(i,j) __mod_vars_MOD_cipw[(i - __mod_vars_MOD_g3l) + j*(__mod_vars_MOD_g3u-__mod_vars_MOD_g3l+1)];

define mpt_ciup(i,j) __mod_vars_MOD_ciup[(i - __mod_vars_MOD_d1l) + j*(__mod_vars_MOD_d1u-__mod_vars_MOD_d1l+1)];
define mpt_civp(i,j) __mod_vars_MOD_civp[(i - __mod_vars_MOD_d2l) + j*(__mod_vars_MOD_d2u-__mod_vars_MOD_d2l+1)];
define mpt_ciwp(i,j) __mod_vars_MOD_ciwp[(i - __mod_vars_MOD_d3l) + j*(__mod_vars_MOD_d3u-__mod_vars_MOD_d3l+1)];
// filter
define mpt_cfp1(i,j) __mod_vars_MOD_cfp1[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cfp2(i,j) __mod_vars_MOD_cfp2[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cfp3(i,j) __mod_vars_MOD_cfp3[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];

define mpt_cfu1(i,j) __mod_vars_MOD_cfu1[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cfv2(i,j) __mod_vars_MOD_cfv2[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cfw3(i,j) __mod_vars_MOD_cfw3[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];
// integrator (not for pressure grid)
define mpt_cint1(i,j) __mod_vars_MOD_cint1[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l+1)];
define mpt_cint2(i,j) __mod_vars_MOD_cint2[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l+1)];
define mpt_cint3(i,j) __mod_vars_MOD_cint3[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l+1)];
// 2nd derivative (multigrid)
define mpt_cp11r(i,j,k) __mod_vars_MOD_cp11r[(i + 1) + j*3 + (k - 1)*3*(__mod_vars_MOD_n1+1)];
define mpt_cp22r(i,j,k) __mod_vars_MOD_cp22r[(i + 1) + j*3 + (k - 1)*3*(__mod_vars_MOD_n2+1)];
define mpt_cp33r(i,j,k) __mod_vars_MOD_cp33r[(i + 1) + j*3 + (k - 1)*3*(__mod_vars_MOD_n3+1)];

define mpt_cu11r(i,j,k) __mod_vars_MOD_cu11r[(i + 1) + j*3 + (k - 1)*3*(__mod_vars_MOD_n1+1)];
define mpt_cv22r(i,j,k) __mod_vars_MOD_cv22r[(i + 1) + j*3 + (k - 1)*3*(__mod_vars_MOD_n2+1)];
define mpt_cw33r(i,j,k) __mod_vars_MOD_cw33r[(i + 1) + j*3 + (k - 1)*3*(__mod_vars_MOD_n3+1)];

define mpt_cdg1(i,j,k) __mod_vars_MOD_cdg1[(i + 1) + (j - 1)*3 + (k - 1)*3*__mod_vars_MOD_n1];
define mpt_cdg2(i,j,k) __mod_vars_MOD_cdg2[(i + 1) + (j - 1)*3 + (k - 1)*3*__mod_vars_MOD_n2];
define mpt_cdg3(i,j,k) __mod_vars_MOD_cdg3[(i + 1) + (j - 1)*3 + (k - 1)*3*__mod_vars_MOD_n3];
// interpolation (multigrid)
define mpt_ci1(i,j,k) __mod_vars_MOD_ci1[(i - 1) + (j - 1)*2 + (k - 1)*2*__mod_vars_MOD_n1];
define mpt_ci2(i,j,k) __mod_vars_MOD_ci2[(i - 1) + (j - 1)*2 + (k - 1)*2*__mod_vars_MOD_n2];
define mpt_ci3(i,j,k) __mod_vars_MOD_ci3[(i - 1) + (j - 1)*2 + (k - 1)*2*__mod_vars_MOD_n3];

define mpt_cih1(i,j) __mod_vars_MOD_cih1[(i - 1) + j*2];
define mpt_cih2(i,j) __mod_vars_MOD_cih2[(i - 1) + j*2];
define mpt_cih3(i,j) __mod_vars_MOD_cih3[(i - 1) + j*2];
// restriction (multigrid)
define mpt_cr1(i,j,k) __mod_vars_MOD_cr1[(i + 1) + (j - 1)*3 + (k - 2)*3*__mod_vars_MOD_n1];
define mpt_cr2(i,j,k) __mod_vars_MOD_cr2[(i + 1) + (j - 1)*3 + (k - 2)*3*__mod_vars_MOD_n2];
define mpt_cr3(i,j,k) __mod_vars_MOD_cr3[(i + 1) + (j - 1)*3 + (k - 2)*3*__mod_vars_MOD_n3];

define mpt_crest1(i,j,k) __mod_vars_MOD_crest1[(i - __mod_vars_MOD_b1l) + j*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - 1)*(__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n1+1)];
define mpt_crest2(i,j,k) __mod_vars_MOD_crest2[(i - __mod_vars_MOD_b2l) + j*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (k - 1)*(__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n2+1)];
define mpt_crest3(i,j,k) __mod_vars_MOD_crest3[(i - __mod_vars_MOD_b3l) + j*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l) + (k - 1)*(__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)*(__mod_vars_MOD_n3+1)];

define mpt_crh1(i,j) __mod_vars_MOD_crh1[(i - 1) + (j - 1)*2];
define mpt_crh2(i,j) __mod_vars_MOD_crh2[(i - 1) + (j - 1)*2];
define mpt_crh3(i,j) __mod_vars_MOD_crh3[(i - 1) + (j - 1)*2];
// grid specifications
// physical coordinates (global)
define mpt_y1p(i) __mod_vars_MOD_y1p[i - 1];
define mpt_y1u(i) __mod_vars_MOD_y1u[i];
	
define mpt_y2p(i) __mod_vars_MOD_y2p[i - 1];
define mpt_y2v(i) __mod_vars_MOD_y2v[i];
	
define mpt_3p(i) __mod_vars_MOD_y3p[i - 1];
define mpt_3w(i) __mod_vars_MOD_y3w[i];
// physical coordinates (block)
define mpt_x1p(i) __mod_vars_MOD_x1p[i - __mod_vars_MOD_b1l];
define mpt_x1u(i) __mod_vars_MOD_x1u[i - __mod_vars_MOD_b1l];

define mpt_x2p(i) __mod_vars_MOD_x2p[i - __mod_vars_MOD_b2l];
define mpt_x2v(i) __mod_vars_MOD_x2v[i - __mod_vars_MOD_b2l];
	
define mpt_x3p(i) __mod_vars_MOD_x3p[i - __mod_vars_MOD_b3l];
define mpt_x3w(i) __mod_vars_MOD_x3w[i - __mod_vars_MOD_b3l];
// physical coordinates (block, multigrid)
define mpt_x1pr(i,j) __mod_vars_MOD_x1pr[(i - __mod_vars_MOD_b1l) + (j - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)];
define mpt_x1ur(i,j) __mod_vars_MOD_x1ur[(i - __mod_vars_MOD_b1l) + (j - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)];
	
define mpt_x2pr(i,j) __mod_vars_MOD_x2pr[(i - __mod_vars_MOD_b2l) + (j - 1)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_x2vr(i,j) __mod_vars_MOD_x2vr[(i - __mod_vars_MOD_b2l) + (j - 1)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
	
define mpt_x3pr(i,j) __mod_vars_MOD_x3pr[(i - __mod_vars_MOD_b3l) + (j - 1)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
define mpt_x3wr(i,j) __mod_vars_MOD_x3wr[(i - __mod_vars_MOD_b3l) + (j - 1)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
// grid widths (global)	
define mpt_dy1p(i) __mod_vars_MOD_dy1p[i - 1];
define mpt_dy1u(i) __mod_vars_MOD_dy1u[i];
	
define mpt_dy2p(i) __mod_vars_MOD_dy2p[i - 1];
define mpt_dy2v(i) __mod_vars_MOD_dy2v[i];
	
define mpt_dy3p(i) __mod_vars_MOD_dy3p[i - 1];
define mpt_dy3w(i) __mod_vars_MOD_dy3w[i];
// grid widths (block)
define mpt_dx1p(i) __mod_vars_MOD_dx1p[i - 1];
define mpt_dx1u(i) __mod_vars_MOD_dx1u[i];

define mpt_dx2p(i) __mod_vars_MOD_dx2p[i - 1];
define mpt_dx2v(i) __mod_vars_MOD_dx2v[i];

define mpt_dx3p(i) __mod_vars_MOD_dx3p[i - 1];
define mpt_dx3w(i) __mod_vars_MOD_dx3w[i];

define mpt_dx1dm(i) __mod_vars_MOD_dx1dm[i - 1];
define mpt_dx1pm(i) __mod_vars_MOD_dx1pm[i - 1];
define mpt_ddx1pm(i) __mod_vars_MOD_ddx1pm[i - 1];

define mpt_dx2dm(i) __mod_vars_MOD_dx2dm[i - 1];
define mpt_dx2pm(i) __mod_vars_MOD_dx2pm[i - 1];
define mpt_ddx2pm(i) __mod_vars_MOD_ddx2pm[i - 1];

define mpt_dx3dm(i) __mod_vars_MOD_dx3dm[i - 1];
define mpt_dx3pm(i) __mod_vars_MOD_dx3pm[i - 1];
define mpt_ddx3pm(i) __mod_vars_MOD_ddx3pm[i - 1];

define mpt_dx1gm(i) __mod_vars_MOD_dx1gm[i];
define mpt_dx1um(i) __mod_vars_MOD_dx1um[i];
define mpt_ddx1um(i) __mod_vars_MOD_ddx1um[i];

define mpt_dx2gm(i) __mod_vars_MOD_dx2gm[i];
define mpt_dx2vm(i) __mod_vars_MOD_dx2vm[i];
define mpt_ddx2vm(i) __mod_vars_MOD_ddx2vm[i];

define mpt_dx3gm(i) __mod_vars_MOD_dx3gm[i];
define mpt_dx3wm(i) __mod_vars_MOD_dx3wm[i];
define mpt_ddx3wm(i) __mod_vars_MOD_ddx3wm[i];
// work fields
// velocities
define mpt_vel(i,j,k,l) __mod_vars_MOD_vel[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (l - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
// non-linear term
define mpt_nl(i,j,k,l) __mod_vars_MOD_nl[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (l - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
// right hand side
define mpt_rhs(i,j,k,l) __mod_vars_MOD_rhs[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (l - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
// pressure
define mpt_pre(i,j,k) __mod_vars_MOD_pre[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
// outflow boundary conditions (velocity field)
define mpt_bc11(i,j,k) __mod_vars_MOD_bc11[(i - 1) + (j - 1)*__mod_vars_MOD_n2 + (k - 1)*__mod_vars_MOD_n2*__mod_vars_mod_n3];
define mpt_nlbc11(i,j,k) __mod_vars_MOD_nlbc11[(i - 1) + (j - 1)*__mod_vars_MOD_n2 + (k - 1)*__mod_vars_MOD_n2*__mod_vars_mod_n3];

define mpt_bc12(i,j,k) __mod_vars_MOD_bc12[i + (j - 1)*(__mod_vars_MOD_n1+1) + (k - 1)*(__mod_vars_MOD_n1+1)*__mod_vars_mod_n3];
define mpt_nlbc12(i,j,k) __mod_vars_MOD_nlbc12[i + (j - 1)*(__mod_vars_MOD_n1+1) + (k - 1)*(__mod_vars_MOD_n1+1)*__mod_vars_mod_n3];

define mpt_bc13(i,j,k) __mod_vars_MOD_bc13[i + (j - 1)*(__mod_vars_MOD_n1+1) + (k - 1)*(__mod_vars_MOD_n1+1)*__mod_vars_mod_n2];
define mpt_nlbc13(i,j,k) __mod_vars_MOD_nlbc13[i + (j - 1)*(__mod_vars_MOD_n1+1) + (k - 1)*(__mod_vars_MOD_n1+1)*__mod_vars_mod_n2];

define mpt_bc21(i,j,k) __mod_vars_MOD_bc21[i + (j - 1)*(__mod_vars_MOD_n2+1) + (k - 1)*(__mod_vars_MOD_n2+1)*__mod_vars_mod_n3];
define mpt_nlbc21(i,j,k) __mod_vars_MOD_nlbc21[i + (j - 1)*(__mod_vars_MOD_n2+1) + (k - 1)*(__mod_vars_MOD_n2+1)*__mod_vars_mod_n3];

define mpt_bc22(i,j,k) __mod_vars_MOD_bc22[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_MOD_n3];
define mpt_nlbc22(i,j,k) __mod_vars_MOD_nlbc22[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_MOD_n3];

define mpt_bc23(i,j,k) __mod_vars_MOD_bc23[(i - 1) + j*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*(__mod_vars_mod_n2+1)];
define mpt_nlbc23(i,j,k) __mod_vars_MOD_nlbc23[(i - 1) + j*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*(__mod_vars_mod_n2+1)];

define mpt_bc31(i,j,k) __mod_vars_MOD_bc31[(i - 1) + j*__mod_vars_MOD_n2 + (k - 1)*__mod_vars_MOD_n2*(__mod_vars_mod_n3+1)];
define mpt_nlbc31(i,j,k) __mod_vars_MOD_nlbc31[(i - 1) + j*__mod_vars_MOD_n2 + (k - 1)*__mod_vars_MOD_n2*(__mod_vars_mod_n3+1)];

define mpt_bc32(i,j,k) __mod_vars_MOD_bc32[(i - 1) + j*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*(__mod_vars_mod_n3+1)];
define mpt_nlbc32(i,j,k) __mod_vars_MOD_nlbc32[(i - 1) + j*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*(__mod_vars_mod_n3+1)];

define mpt_bc33(i,j,k) __mod_vars_MOD_bc33[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_mod_n2];
define mpt_nlbc33(i,j,k) __mod_vars_MOD_nlbc33[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_mod_n2];

define mpt_drift1(i,j,k) __mod_vars_MOD_drift1[(i - __mod_vars_MOD_b2l) + (j - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (k - 1)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

define mpt_drift2(i,j,k) __mod_vars_MOD_drift2[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

define mpt_drift3(i,j,k) __mod_vars_MOD_drift3[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

// residual
define mpt_res(i,j,k) __mod_vars_MOD_res[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
// pressure gradient
define mpt_gpre(i,j,k) __mod_vars_MOD_gpre[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
// weights for divergence-freeness
define mpt_weight(i,j,k) __mod_vars_MOD_weight[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_MOD_n2];
// null-space vector	
define mpt_psi __mod_vars_MOD_psi[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_psi_vel __mod_vars_MOD_psi_vel[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (l - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];

define mpt_psi_rel1 __mod_vars_MOD_psi_rel1[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
	// NOTE: for the variables below nn is needed and has not yet been interfaced (may need to be moved above these lines)
define mpt_psi_rel2(i,j,k) __mod_vars_MOD_psi_rel2[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel3(i,j,k) __mod_vars_MOD_psi_rel3[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel4(i,j,k) __mod_vars_MOD_psi_rel4[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel5(i,j,k) __mod_vars_MOD_psi_rel5[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel6(i,j,k) __mod_vars_MOD_psi_rel6[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel7(i,j,k) __mod_vars_MOD_psi_rel7[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel8(i,j,k) __mod_vars_MOD_psi_rel8[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel9(i,j,k) __mod_vars_MOD_psi_rel9[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel10(i,j,k) __mod_vars_MOD_psi_rel10[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel11(i,j,k) __mod_vars_MOD_psi_rel11[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel12(i,j,k) __mod_vars_MOD_psi_rel12[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel13(i,j,k) __mod_vars_MOD_psi_rel13[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel14(i,j,k) __mod_vars_MOD_psi_rel14[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_psi_rel15(i,j,k) __mod_vars_MOD_psi_rel15[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*14]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
	
define mpt_th11(i,j,k) __mod_vars_MOD_th11[(i - 1) + (j - 1)*__mod_vars_MOD_n2 + (k - 1)*__mod_vars_MOD_n2*__mod_vars_mod_n3];

define mpt_th12(i,j,k) __mod_vars_MOD_th12[i + (j - 1)*(__mod_vars_MOD_n1+1) + (k - 1)*(__mod_vars_MOD_n1+1)*__mod_vars_mod_n3];

define mpt_th13(i,j,k) __mod_vars_MOD_th13[i + (j - 1)*(__mod_vars_MOD_n1+1) + (k - 1)*(__mod_vars_MOD_n1+1)*__mod_vars_mod_n2];

define mpt_th21(i,j,k) __mod_vars_MOD_th21[i + (j - 1)*(__mod_vars_MOD_n2+1) + (k - 1)*(__mod_vars_MOD_n2+1)*__mod_vars_mod_n3];

define mpt_th22(i,j,k) __mod_vars_MOD_th22[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_MOD_n3];

define mpt_th23(i,j,k) __mod_vars_MOD_th23[(i - 1) + j*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*(__mod_vars_mod_n2+1)];

define mpt_th31(i,j,k) __mod_vars_MOD_th31[(i - 1) + j*__mod_vars_MOD_n2 + (k - 1)*__mod_vars_MOD_n2*(__mod_vars_mod_n3+1)];

define mpt_th32(i,j,k) __mod_vars_MOD_th32[(i - 1) + j*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*(__mod_vars_mod_n3+1)];

define mpt_th33(i,j,k) __mod_vars_MOD_th33[(i - 1) + (j - 1)*__mod_vars_MOD_n1 + (k - 1)*__mod_vars_MOD_n1*__mod_vars_mod_n2];


// multigrid
define mpt_vec1c(i,j,k) __mod_vars_MOD_vec1c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec2a(i,j,k) __mod_vars_MOD_vec2a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec2b(i,j,k) __mod_vars_MOD_vec2b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec2c(i,j,k) __mod_vars_MOD_vec2c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*1]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*1]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec3a(i,j,k) __mod_vars_MOD_vec3a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec3b(i,j,k) __mod_vars_MOD_vec3b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec3b(i,j,k) __mod_vars_MOD_vec3c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*2]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*2]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec4a(i,j,k) __mod_vars_MOD_vec4a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec4b(i,j,k) __mod_vars_MOD_vec4b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec4c(i,j,k) __mod_vars_MOD_vec4c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*3]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*3]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec5a(i,j,k) __mod_vars_MOD_vec5a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec5b(i,j,k) __mod_vars_MOD_vec5b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec5c(i,j,k) __mod_vars_MOD_vec5c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*4]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*4]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec6a(i,j,k) __mod_vars_MOD_vec6a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec6b(i,j,k) __mod_vars_MOD_vec6b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec6c(i,j,k) __mod_vars_MOD_vec6c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*5]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*5]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec7a(i,j,k) __mod_vars_MOD_vec7a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec7b(i,j,k) __mod_vars_MOD_vec7b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec7c(i,j,k) __mod_vars_MOD_vec7c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*6]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*6]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec8a(i,j,k) __mod_vars_MOD_vec8a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec8b(i,j,k) __mod_vars_MOD_vec8b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec8c(i,j,k) __mod_vars_MOD_vec8c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*7]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*7]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec9a(i,j,k) __mod_vars_MOD_vec9a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec9b(i,j,k) __mod_vars_MOD_vec9b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec9c(i,j,k) __mod_vars_MOD_vec9c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*8]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*8]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec10a(i,j,k) __mod_vars_MOD_vec10a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec10b(i,j,k) __mod_vars_MOD_vec10b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec10c(i,j,k) __mod_vars_MOD_vec10c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*9]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*9]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec11a(i,j,k) __mod_vars_MOD_vec11a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec11b(i,j,k) __mod_vars_MOD_vec11b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec11c(i,j,k) __mod_vars_MOD_vec11c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*10]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*10]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec12a(i,j,k) __mod_vars_MOD_vec12a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec12b(i,j,k) __mod_vars_MOD_vec12b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec12c(i,j,k) __mod_vars_MOD_vec12c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*11]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*11]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec13a(i,j,k) __mod_vars_MOD_vec13a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec13b(i,j,k) __mod_vars_MOD_vec13b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec13c(i,j,k) __mod_vars_MOD_vec13c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*12]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*12]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec14a(i,j,k) __mod_vars_MOD_vec14a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec14b(i,j,k) __mod_vars_MOD_vec14b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec14c(i,j,k) __mod_vars_MOD_vec14c[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*13]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*13]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_vec15a(i,j,k) __mod_vars_MOD_vec15a[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*14]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
define mpt_vec15b(i,j,k) __mod_vars_MOD_vec15b[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_nn[0+3*14]+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_nn[1+3*14]+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

// BiCGstab / Richardson
define mpt_pp(i,j,k) __mod_vars_MOD_pp[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_ap(i,j,k) __mod_vars_MOD_ap[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_rr(i,j,k) __mod_vars_MOD_rr[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_rh(i,j,k) __mod_vars_MOD_rh[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_ar(i,j,k) __mod_vars_MOD_ar[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_z1(i,j,k) __mod_vars_MOD_z1[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_z2(i,j,k) __mod_vars_MOD_z2[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
// product_div_grad
define mpt_dig(i,j,k) __mod_vars_MOD_dig[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
// helper fields (pressure iteration)
define mpt_work1(i,j,k) __mod_vars_MOD_work1[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_work2(i,j,k) __mod_vars_MOD_work2[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];

define mpt_work3(i,j,k) __mod_vars_MOD_work3[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)];
// line relaxation
define mpt_vec1(i) __mod_vars_MOD_vec1[i - 1];
define mpt_dia1(i) __mod_vars_MOD_dia1[i - 1];
define mpt_sor1(i) __mod_vars_MOD_sor1[i - 1];
define mpt_band1(i,j) __mod_vars_MOD_band1[(i - 1) + (j - 1)*2];

define mpt_vec2(i) __mod_vars_MOD_vec2[i - 1];
define mpt_dia2(i) __mod_vars_MOD_dia2[i - 1];
define mpt_sor2(i) __mod_vars_MOD_sor2[i - 1];
define mpt_band2(i,j) __mod_vars_MOD_band2[(i - 1) + (j - 1)*2];

define mpt_vec3(i) __mod_vars_MOD_vec3[i - 1];
define mpt_dia3(i) __mod_vars_MOD_dia3[i - 1];
define mpt_sor3(i) __mod_vars_MOD_sor3[i - 1];
define mpt_band3(i,j) __mod_vars_MOD_band3[(i - 1) + (j - 1)*2];

// indexing (interval limits, shifts)
// index shift (block -> global)
define mpt_ishift __mod_vars_MOD_ishift;
define mpt_jshift __mod_vars_MOD_jshift;
define mpt_kshift __mod_vars_MOD_kshift;
// domain size (cleaned from periodicity)
define mpt_dim1 __mod_vars_MOD_dim1;
define mpt_dim2 __mod_vars_MOD_dim2;
define mpt_dim3 __mod_vars_MOD_dim3;
// pressure / concentrations (incl. border)
define mpt_s1p __mod_vars_MOD_s1p;
define mpt_s2p __mod_vars_MOD_s2p;
define mpt_s3p __mod_vars_MOD_s3p;

define mpt_n1p __mod_vars_MOD_n1p;
define mpt_n2p __mod_vars_MOD_n2p;
define mpt_n3p __mod_vars_MOD_n3p;
// velocities (incl. border)
define mpt_s11b __mod_vars_MOD_s11b;
define mpt_s21b __mod_vars_MOD_s21b;
define mpt_s31b __mod_vars_MOD_s31b;
	
define mpt_s12b __mod_vars_MOD_s12b;
define mpt_s22b __mod_vars_MOD_s22b;
define mpt_s32b __mod_vars_MOD_s32b;
	
define mpt_s13b __mod_vars_MOD_s13b;
define mpt_s23b __mod_vars_MOD_s23b;
define mpt_s33b __mod_vars_MOD_s33b;

define mpt_n11b __mod_vars_MOD_n11b;
define mpt_n21b __mod_vars_MOD_n21b;
define mpt_n31b __mod_vars_MOD_n31b;
	
define mpt_n12b __mod_vars_MOD_n12b;
define mpt_n22b __mod_vars_MOD_n22b;
define mpt_n32b __mod_vars_MOD_n32b;
	
define mpt_n13b __mod_vars_MOD_n13b;
define mpt_n23b __mod_vars_MOD_n23b;
define mpt_n33b __mod_vars_MOD_n33b;
// velocities (excl. border)
define mpt_s11 __mod_vars_MOD_s11;
define mpt_s21 __mod_vars_MOD_s21;
define mpt_s31 __mod_vars_MOD_s31;
	
define mpt_s12 __mod_vars_MOD_s12;
define mpt_s22 __mod_vars_MOD_s22;
define mpt_s32 __mod_vars_MOD_s32;
	
define mpt_s13 __mod_vars_MOD_s13;
define mpt_s23 __mod_vars_MOD_s23;
define mpt_s33 __mod_vars_MOD_s33;

define mpt_n11 __mod_vars_MOD_n11;
define mpt_n21 __mod_vars_MOD_n21;
define mpt_n31 __mod_vars_MOD_n31;
	
define mpt_n12 __mod_vars_MOD_n12;
define mpt_n22 __mod_vars_MOD_n22;
define mpt_n32 __mod_vars_MOD_n32;
	
define mpt_n13 __mod_vars_MOD_n13;
define mpt_n23 __mod_vars_MOD_n23;
define mpt_n33 __mod_vars_MOD_n33;
// coarse grids (multigrid, incl. border)
define mpt_s1r __mod_vars_MOD_s1r;
define mpt_s2r __mod_vars_MOD_s2r;
define mpt_s3r __mod_vars_MOD_s3r;

define mpt_d1r __mod_vars_MOD_d1r;
define mpt_d2r __mod_vars_MOD_d2r;
define mpt_d3r __mod_vars_MOD_d3r;
// coarse grids (multigrid, excl. border)
define mpt_s11r __mod_vars_MOD_s11r;
define mpt_s22r __mod_vars_MOD_s22r;
define mpt_s33r __mod_vars_MOD_s33r;

define mpt_d11r __mod_vars_MOD_d11r;
define mpt_d22r __mod_vars_MOD_d22r;
define mpt_d33r __mod_vars_MOD_d33r;
// overlapping convention of the blocks (multigrid, see mod_setup)
define mpt_ls1 __mod_vars_MOD_ls1;
define mpt_ls2 __mod_vars_MOD_ls2;
define mpt_ls3 __mod_vars_MOD_ls3;
// exchange direction (multigrid)
define mpt_ex1 __mod_vars_MOD_ex1;
define mpt_ex2 __mod_vars_MOD_ex2;
define mpt_ex3 __mod_vars_MOD_ex3;

// boundary conditions
// global
define mpt_outlet(i,j,k) __mod_vars_MOD_outlet[(i - 1) + (j - 1)*3 + (k - 1)*3*2];

define mpt_bc_1l_global __mod_vars_MOD_bc_1l_global;
define mpt_bc_1u_global __mod_vars_MOD_bc_1u_global;
define mpt_bc_2l_global __mod_vars_MOD_bc_2l_global;
define mpt_bc_2u_global __mod_vars_MOD_bc_2u_global;
define mpt_bc_3l_global __mod_vars_MOD_bc_3l_global;
define mpt_bc_3u_global __mod_vars_MOD_bc_3u_global;
// local (block)
define mpt_bc_1l __mod_vars_MOD_bc_1l;
define mpt_bc_1u __mod_vars_MOD_bc_1u;
define mpt_bc_2l __mod_vars_MOD_bc_2l;
define mpt_bc_2u __mod_vars_MOD_bc_2u;
define mpt_bc_3l __mod_vars_MOD_bc_3l;
define mpt_bc_3u __mod_vars_MOD_bc_3u;
// field properties
define mpt_n_gather(i,j) __mod_vars_MOD_n_gather[(i - 1) + (j - 1)*3];
define mpt_nn(i,j) __mod_vars_MOD_nn[(i - 1) + (j - 1)*3];
define mpt_nb(i,j) __mod_vars_MOD_nb[(i - 1) + (j - 1)*3];
define mpt_ib(i,j) __mod_vars_MOD_ib[(i - 1) + (j - 1)*3];
define mpt_snf(i,j,k) __mod_vars_MOD_snf[(i - 1) + (j - 1)*2 + (k - 1)*2*3];
define mpt_snb(i,j,k) __mod_vars_MOD_snb[(i - 1) + (j - 1)*2 + (k - 1)*2*3];
define mpt_bc(i,j,k) __mod_vars_MOD_bc[(i - 1) + (j - 1)*2 + (k - 1)*2*3];
define mpt_ngb(i,j,k) __mod_vars_MOD_ngb[(i - 1) + (j - 1)*2 + (k - 1)*2*3];
define mpt_comm1(i) __mod_vars_MOD_comm1[i - 1];
define mpt_comm2(i) __mod_vars_MOD_comm2[i - 1];
define mpt_rankc2(i) __mod_vars_MOD_rankc2[i - 1];
define mpt_participate_yes(i) __mod_vars_MOD_participate_yes[i - 1];
define mpt_recvr(i,j) __mod_vars_MOD_recvr[(i - 1) + (j - 1)*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_recvi(i,j) __mod_vars_MOD_recvi[(i - 1) + (j - 1)*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_dispr(i,j) __mod_vars_MOD_dispr[(i - 1) + (j - 1)*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_dispi(i,j) __mod_vars_MOD_dispi[(i - 1) + (j - 1)*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_offsr(i,j,k) __mod_vars_MOD_offsr[(i - 1) + (j - 1)*3 + (k - 1)*3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_offsi(i,j,k) __mod_vars_MOD_offsi[(i - 1) + (j - 1)*3 + (k - 1)*3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_sizsr(i,j,k) __mod_vars_MOD_sizsr[(i - 1) + (j - 1)*3 + (k - 1)*3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
define mpt_sizsi(i,j,k) __mod_vars_MOD_sizsi[(i - 1) + (j - 1)*3 + (k - 1)*3*__mod_dims_MOD_nb1*__mod_dims_MOD_nb2*__mod_dims_MOD_nb3];
// physical parameters
define mpt_l1 __mod_vars_MOD_l1;
define mpt_l2 __mod_vars_MOD_l2;
define mpt_l3 __mod_vars_MOD_l3;
define mpt_re __mod_vars_MOD_re;
// numerical parameters
// general
define mpt_cfl __mod_vars_MOD_cfl;
define mpt_time __mod_vars_MOD_time;
define mpt_dtime __mod_vars_MOD_dtime;
define mpt_subtime __mod_vars_MOD_subtime;
define mpt_time_start __mod_vars_MOD_time_start;
define mpt_time_end __mod_vars_MOD_time_end;
define mpt_dtime_max __mod_vars_MOD_dtime_max;
define mpt_dtime0 __mod_vars_MOD_dtime0;
define mpt_dtime_old __mod_vars_MOD_dtime_old;
define mpt_timestep __mod_vars_MOD_timestep;
define mpt_timestep_old __mod_vars_MOD_timestep_old;
define mpt_substep __mod_vars_MOD_substep;
define mpt_n_timesteps __mod_vars_MOD_n_timesteps;
define mpt_mapping_yes __mod_vars_MOD_mapping_yes;
define mpt_upwind_yes __mod_vars_MOD_upwind_yes;
define mpt_euler_yes __mod_vars_MOD_euler_yes;
define mpt_stokes_yes __mod_vars_MOD_stokes_yes;
define mpt_twpstep_yes __mod_vars_MOD_twostep_yes;
define mpt_filter_bc_yes __mod_vars_MOD_filter_bc_yes;
define mpt_timeint_mode __mod_vars_MOD_timeint_mode;
define mpt_forcing_mode __mod_vars_MOD_forcing_mode;
define mpt_bulkflow_dir __mod_vars_MOD_bulkflow_dir;
define mpt_n_lp_vel __mod_vars_MOD_n_lp_vel;
define mpt_n_hp_vel __mod_vars_MOD_n_hp_vel;
define mpt_chi_vel __mod_vars_MOD_chi_vel;
define mpt_vel_bulk __mod_vars_MOD_vel_bulk;
// Runge-Kutta coefficients
define mpt_ark(i) __mod_vars_MOD_ark[i - 1];
define mpt_brk(i) __mod_vars_MOD_brk[i - 1];
define mpt_rk_steps __mod_vars_MOD_rk_steps;
// look-up table for stability region of time-integration
define mpt_stabilitylimit(i) __mod_vars_MOD_stabilitylimit[i];
define mpt_n_stab __mod_vars_MOD_n_stab;
// Helmholtz pre-factors
define mpt_thetal __mod_vars_MOD_thetal;
define mpt_multl __mod_vars_MOD_multl;
// temporal control
define mpt_int_dtime __mod_vars_MOD_int_dtime;
define mpt_int_lev_pre __mod_vars_MOD_int_lev_pre;
	
define mpt_stride_large __mod_vars_MOD_stride_large[i - 1];
define mpt_stride_med __mod_vars_MOD_stride_med[i - 1];
define mpt_stride_small __mod_vars_MOD_stride_small[i - 1];
define mpt_write_large __mod_vars_MOD_write_large;
define mpt_write_med __mod_vars_MOD_write_med;
define mpt_write_small __mod_vars_MOD_write_small;
define mpt_time_out_scal __mod_vars_MOD_time_out_scal;
define mpt_dtime_out_scal __mod_vars_MOD_dtime_out_scal;
define mpt_time_out_vect __mod_vars_MOD_time_out_vect;
define mpt_dtime_out_vect __mod_vars_MOD_dtime_out_vect;
	
define mpt_write_out_scal __mod_vars_MOD_write_out_scal;
define mpt_write_out_vect __mod_vars_MOD_write_out_vect;
define mpt_new_dtime __mod_vars_MOD_new_dtime;
define mpt_finish_yes __mod_vars_MOD_finish_yes;

define mpt_write_count __mod_vars_MOD_write_count;
define mpt_restart __mod_vars_MOD_restart;
define mpt_restart_char __mod_vars_MOD_restart_char; // has length 3
define mpt_n_conc_old __mod_vars_MOD_n_conc_old;
// further control options
define mpt_task __mod_vars_MOD_task;
define mpt_read_nullspace_yes __mod_vars_MOD_read_nullspace_yes;
define mpt_nullspace_yes __mod_vars_MOD_nullspace_yes;
define mpt_nullspace_coarse_yes __mod_vars_MOD_nullspace_coarse_yes;
define mpt_nullspace_ortho_yes __mod_vars_MOD_nullspace_ortho_yes;

define mpt_write_stout_yes __mod_vars_MOD_write_stout_yes;
define mpt_log_iteration_yes __mod_vars_MOD_log_iteration_yes;
define mpt_write_restart_yes __mod_vars_MOD_write_restart_yes;
define mpt_write_lambda2_yes __mod_vars_MOD_write_lambda2_yes;
define mpt_write_test_yes __mod_vars_MOD_write_test_yes;
// global running indices
define mpt_direction __mod_vars_MOD_direction;
// explicit treatment of corners with Dirichlet boundary conditions
define mpt_corner_yes __mod_vars_MOD_corner_yes;
// system time
define mpt_elatime __mod_vars_MOD_elatime;
define mpt_ctime(i)__mod_vars_MOD_ctime[i - 1];
define mpt_day __mod_vars_MOD_day;
define mpt_hour __mod_vars_MOD_hour;
define mpt_minu __mod_vars_MOD_minu;
define mpt_sec __mod_vars_MOD_sec;
define mpt_msec __mod_vars_MOD_msec;
// iteration parameters
// termination criterion / absolute accuracy of velocities
define mpt_epsu __mod_vars_MOD_epsu;
define mpt_epsu0 __mod_vars_MOD_epsu0;
// smoother
define mpt_jacobi_yes __mod_vars_MOD_jacobi_yes;
// number of max. iterations
define mpt_n_it_outer __mod_vars_MOD_n_it_outer;
define mpt_n_it_poisson __mod_vars_MOD_n_it_poisson;
define mpt_n_it_helmh_vel __mod_vars_MOD_n_it_helmh_vel;
// expected convergence rate (outer iteration)
define mpt_precratio0(i) __mod_vars_MOD_precratio0[i - 1];
define mpt_precoffset(i) __mod_vars_MOD_precoffset0[i - 1];
define mpt_precratio(i,j) __mod_vars_MOD_precratio[(i - 1) + (j - 1)*__mod_vars_MOD_rk_steps];
define mpt_precoffset(i,j) __mod_vars_MOD_precoffset[(i - 1) + (j - 1)*__mod_vars_MOD_rk_steps];
// null-initialization (outer iteration)
define mpt_init_pre __mod_vars_MOD_init_pre[i - 1];
define mpt_init_vel __mod_vars_MOD_init_vel[i - 1];
// preconditioning (multigrid)
define mpt_precond_outer __mod_vars_MOD_precond_outer;
define mpt_precond_poisson __mod_vars_MOD_precond_poisson;
define mpt_precond_helmh_vel __mod_vars_MOD_precond_helmh_vel;
// number of smoothings per grid level (multigrid)
define mpt_n_relax_down __mod_vars_MOD_n_relax_down;
define mpt_n_relax_up __mod_vars_MOD_n_relax_up;
define mpt_n_relax_bottom __mod_vars_MOD_n_relax_bottom;
// implicit directions of line relaxations (multigrid)
define mpt_impl_dir(i) __mod_vars_MOD_impl_dir[i - 1];
	
define mpt_weighting_yes __mod_vars_MOD_weighting_yes;
	
// iteration stats
define mpt_dtime_average __mod_vars_MOD_dtime_average;
define mpt_max_div_init(i) __mod_vars_MOD_max_div_init[i - 1];
define mpt_number_poisson __mod_vars_MOD_number_poisson;
// counters
define mpt_count0(i) __mod_vars_MOD_count0[i - 1];
define mpt_countp(i,j) __mod_vars_MOD_countp[(i - 1) + (j - 1)*__mod_vars_MOD_rk_steps];
define mpt_counth(i,j) __mod_vars_MOD_counth[(i - 1) + (j - 1)*__mod_vars_MOD_rk_steps];
// convergence rate
define mpt_ratio(i) __mod_vars_MOD_ratio0[i - 1];
define mpt_ratiop(i,j) __mod_vars_MOD_ratiop[(i - 1) + (j - 1)*__mod_vars_MOD_rk_steps];
define mpt_ratioh(i,j) __mod_vars_MOD_ratioh[(i - 1) + (j - 1)*__mod_vars_MOD_rk_steps];
	
// MPI
// communicators
define mpt_comm_cart __mod_vars_MOD_comm_cart;
	
define mpt_comm_slice1 __mod_vars_MOD_comm_slice1;
define mpt_comm_bar1 __mod_vars_MOD_comm_bar1;
define mpt_comm_slice2 __mod_vars_MOD_comm_slice2;
define mpt_comm_bar2 __mod_vars_MOD_comm_bar2;
define mpt_comm_slice3 __mod_vars_MOD_comm_slice3;
define mpt_comm_bar3 __mod_vars_MOD_comm_bar3;
// dimension and position of the blocks inside the communicators (grid indices)
define mpt_bar1_size(i) __mod_vars_MOD_bar1_size[i - 1];
define mpt_bar1_offset(i) __mod_vars_MOD_bar1_offset[i - 1];
define mpt_bar2_size(i) __mod_vars_MOD_bar2_size[i - 1];
define mpt_bar2_offset(i) __mod_vars_MOD_bar2_offset[i - 1];
define mpt_bar3_size(i) __mod_vars_MOD_bar3_size[i - 1];
define mpt_bar3_offset(i) __mod_vars_MOD_bar3_offset[i - 1];
// ranks of processes
define mpt_rank __mod_vars_MOD_rank;
define mpt_rank_bar1 __mod_vars_MOD_rank_bar1;
define mpt_rank_slice1 __mod_vars_MOD_rank_slice1;
define mpt_rank_bar2 __mod_vars_MOD_rank_bar2;
define mpt_rank_slice2 __mod_vars_MOD_rank_slice2;
define mpt_rank_bar3 __mod_vars_MOD_rank_bar3;
define mpt_rank_slice3 __mod_vars_MOD_rank_slice3;
// ranks of neighboring processes (in cartesian grid)
define mpt_rank1l __mod_vars_MOD_rank1l;
define mpt_rank1u __mod_vars_MOD_rank1u;
define mpt_rank2l __mod_vars_MOD_rank2l;
define mpt_rank2u __mod_vars_MOD_rank2u;
define mpt_rank3l __mod_vars_MOD_rank3l;
define mpt_rank3u __mod_vars_MOD_rank3u;
// error handle
define mpt_merror __mod_vars_MOD_merror;
// request handles
define mpt_req1l __mod_vars_MOD_req1l;
define mpt_req1u __mod_vars_MOD_req1u;
define mpt_req2l __mod_vars_MOD_req2l;
define mpt_req2u __mod_vars_MOD_req2u;
define mpt_req3l __mod_vars_MOD_req3l;
define mpt_req3u __mod_vars_MOD_req3u;
// HDF5
define mpt_herror __mod_vars_MOD_herror;
//*** bbecsek 2015: IMMERSED BOUNDARY PARAMETERS ***
define mpt_reach __mod_vars_MOD_reach;
define mpt_m_bound __mod_vars_MOD_m_bound;
define mpt_m_elems __mod_vars_MOD_m_elems;
define mpt_yb(i,j) __mod_vars_MOD_yb[(i - 1) + (j - 1)*__mod_vars_MOD_m_bound];
define mpt_xb(i,j) __mod_vars_MOD_xb[(i - 1) + (j - 1)*__mod_vars_MOD_m_bound];
define mpt_db(i,j) __mod_vars_MOD_db[(i - 1) + (j - 1)*__mod_vars_MOD_m_bound];
define mpt_ub(i,j) __mod_vars_MOD_ub[(i - 1) + (j - 1)*__mod_vars_MOD_m_bound];
define mpt_fb(i,j) __mod_vars_MOD_fb[(i - 1) + (j - 1)*__mod_vars_MOD_m_bound];
define mpt_container(i,j) __mod_vars_MOD_container[(i - 1) + (j - 1)*__mod_vars_MOD_m_bound];
define mpt_fd(i,j,k,l) __mod_vars_MOD_fd[(i - __mod_vars_MOD_b1l) + (j - __mod_vars_MOD_b2l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l) + (k - __mod_vars_MOD_b3l)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l) + (l - 1)*(__mod_vars_MOD_n1+__mod_vars_MOD_b1u-__mod_vars_MOD_b1l)*(__mod_vars_MOD_n2+__mod_vars_MOD_b2u-__mod_vars_MOD_b2l)*(__mod_vars_MOD_n3+__mod_vars_MOD_b3u-__mod_vars_MOD_b3l)];
 define mpt_ddf_type __mod_vars_MOD_ddf_type;
define mpt_ib_on __mod_vars_MOD_ib_on;
define mpt_mu_fluid __mod_vars_MOD_mu_fluid;
define mpt_rho_fluid __mod_vars_MOD_rho_fluid;
define mpt_e_mod __mod_vars_MOD_e_mod;
define mpt_nu_poiss __mod_vars_MOD_nu_poiss;
define mpt_rho_solid __mod_vars_MOD_rho_solid;
define mpt_strain(i,j) __mod_vars_MOD_strain[(i - 1) + (j - 1)*__mod_vars_MOD_m_elems];
define mpt_stress(i,j) __mod_vars_MOD_stress[(i - 1) + (j - 1)*__mod_vars_MOD_m_elems];
define mpt_k_stiff(i,j) __mod_vars_MOD_k_stiff[(i - 1) + (j - 1)*(2*__mod_vars_MOD_m_bound)];
define mpt_m_mass(i,j) __mod_vars_MOD_m_mass[(i - 1) + (j - 1)*(2*__mod_vars_MOD_m_bound)];
define mpt_c_damp(i,j) __mod_vars_MOD_c_damp[(i - 1) + (j - 1)*(2*__mod_vars_MOD_m_bound)];
define mpt_c(i,j) __mod_vars_MOD_c[(i - 1) + (j - 1)*3];
define mpt_elems(i,j) __mod_vars_MOD_elems[(i - 1) + (j - 1)*__mod_vars_MOD_m_elems];
define mpt_l_ref __mod_vars_MOD_l_ref;
define mpt_u_ref __mod_vars_MOD_u_ref;
define mpt_req_bcast(i) __mod_vars_MOD_req_bcast[i - 1];
define mpt_req_red(i) __mod_vars_MOD_req_red[i - 1];
define mpt_excl_p(i) __mod_vars_MOD_excl_p[i - 1];
define mpt_delta_b1 __mod_vars_MOD_delta_b1;
define mpt_delta_b2 __mod_vars_MOD_delta_b2;
define mpt_delta_b3 __mod_vars_MOD_delta_b3;
define mpt_total_ranks __mod_vars_MOD_total_ranks;
define mpt_m_ebcs __mod_vars_MOD_m_ebcs;
define mpt_ebcs(i,j) __mod_vars_MOD_ebcs[(i - 1) + (j - 1)*__mod_vars_MOD_m_ebcs];
define mpt_node_vol(i) __mod_vars_MOD_node_vol[i - 1];
define mpt_comm_local __mod_vars_MOD_comm_local;
define mpt_comm_inter __mod_vars_MOD_comm_inter;
define mpt_max_ev_freq __mod_vars_MOD_max_ev_freq;
