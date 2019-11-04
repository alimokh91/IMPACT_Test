#ifndef MOD_INOUT_WRAPPER_H
#define MOD_INOUT_WRAPPER_H
#include<cstring>
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void __mod_inout_MOD_write_fields(void);
       void __mod_inout_MOD_write_restart(void);
       void __mod_inout_MOD_read_restart(void);
       void __mod_inout_MOD_write_hdf(char* filename, char* dsetname, 
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3,
		       int* vel_dir, int* stride, double* phi, 
		       int filename_length, int dsetname_length);
       void __mod_inout_MOD_read_hdf(char* filename, char* dsetname, 
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3,
		       int* vel_dir, double* phi, 
		       int filename_length, int dsetname_length); 
       void __mod_inout_MOD_read2_hdf(char* filename, char* dsetname, 
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3,
		       int* vel_dir, double* phi, 
		       int filename_length, int dsetname_length);
       void __mod_inout_MOD_write_hdf_velall(char* filename, 
		       int* ss1, int* ss2, int* ss3, 
		       int* nn1, int* nn2, int* nn3,
		       int* vel_dir, double* phi, 
		       int filename_length);  
       void __mod_inout_MOD_write_2d_hdf(char* filename, char* dsetname,
		       int* n1, int* n2, 
		       int* ss1, int* ss2, 
		       int* nn1, int* nn2,
		       int* ishift, int* jshift, 
		       int* dir, double* phi, 
		       int filename_length, int dsetname_length);
       void __mod_inout_MOD_read2_2d_hdf(char* filename, char* dsetname, 
		       int* n1, int* n2, 
		       int* ss1, int* ss2, 
		       int* nn1, int* nn2,
		       int* ishift, int* jshift, 
		       int* dir, double* phi, 
		       int filename_length, int dsetname_length);
       void __mod_inout_MOD_write_1d_hdf(char* filename, char* dsetname, 
		       int* nn, int* ssp, int* nnp, int* ishift,
		       int* dir, double* phi, 
		       int filename_length, int dsetname_length);
       void __mod_inout_MOD_read2_1d_hdf(char* filename, char* dsetname, 
		       int* nn, int* ssp, int* nnp, int* ishift,
		       int* dir, double* phi, 
		       int filename_length, int dsetname_length);
       void __mod_inout_MOD_write_stats_hdf_4d(char* filename, char* dsetname,
		       int* ss1, int* ss2, int* ss3, int* ss4, 
		       int* nn1, int* nn2, int* nn3, int* nn4, 
		       double* phi, int filename_length, int dsetname_length);
       void __mod_inout_MOD_read_stats_hdf_4d(char* filename, char* dsetname, 
		       int* ss1, int* ss2, int* ss3, int* ss4, 
		       int* nn1, int* nn2, int* nn3, int* nn4, 
		       double* phi, int filename_length, int dsetname_length);             
       void __mod_inout_MOD_start_hdf5_for_testing(void);
       void __mod_inout_MOD_stop_hdf5_for_testing(void);
       void __mod_inout_MOD_filespace_props(int* vel_dir, int* stride, 
		       int* s1w, int* s2w, int* s3w, 
		       int* m1w, int* m2w, int* m3w,
		       int* dim1, int* dim2, int* dim3);
       void __mod_inout_MOD_lambda(void);
       void __mod_inout_MOD_lambda_2d(void);
       void __mod_inout_MOD_write_xdmf_xml(void);
       void __mod_inout_MOD_write_xdmf_timecollection(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_write_fields(){__mod_inout_MOD_write_fields();}
inline void mpt_write_restart(){__mod_inout_MOD_write_restart();}
inline void mpt_read_restart(){__mod_inout_MOD_read_restart();}
inline void mpt_write_hdf(char* filename, char* dsetname, 
		int& ss1, int& ss2, int& ss3,
		int& nn1, int& nn2, int& nn3,
		int& vel_dir, int* stride, double* phi){
	__mod_inout_MOD_write_hdf(filename, dsetname, 
			&ss1, &ss2, &ss3,
			&nn1, &nn2, &nn3, 
			&vel_dir, stride, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_read_hdf(char* filename, char* dsetname, 
		int& ss1, int& ss2, int& ss3, 
		int& nn1, int& nn2, int& nn3,
		int& vel_dir, double* phi){
	__mod_inout_MOD_read_hdf(filename, dsetname, 
			&ss1, &ss2, &ss3,
			&nn1, &nn2, &nn3, 
			&vel_dir, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_read2_hdf(char* filename, char* dsetname, 
		int& ss1, int& ss2, int& ss3, 
		int& nn1, int& nn2, int& nn3,
		int& vel_dir, double* phi){
	__mod_inout_MOD_read2_hdf(filename, dsetname, 
			&ss1, &ss2, &ss3,
			&nn1, &nn2, &nn3, 
			&vel_dir, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_write_hdf_velall(char* filename, 
		int& ss1, int& ss2, int& ss3, 
		int& nn1, int& nn2, int& nn3,
		int& vel_dir, double* phi){
	__mod_inout_MOD_write_hdf_velall(filename, 
			&ss1, &ss2, &ss3,
			&nn1, &nn2, &nn3, 
			&vel_dir, phi, strlen(filename));
}
inline void mpt_write_2d_hdf(char* filename, char* dsetname, 
		int& n1, int& n2, int& ss1, int& ss2, 
		int& nn1, int& nn2,
		int& ishift, int& jshift, 
		int& dir, double* phi){
	__mod_inout_MOD_write_2d_hdf(filename, dsetname, 
			&n1, &n2, &ss1, &ss2, 
			&nn1, &nn2, 
			&ishift, &jshift, 
			&dir, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_read2_2d_hdf(char* filename, char* dsetname, 
		int& n1, int& n2, int& ss1, int& ss2, 
		int& nn1, int& nn2,
		int& ishift, int& jshift, 
		int& dir, double* phi){
	__mod_inout_MOD_read2_2d_hdf(filename, dsetname, 
			&n1, &n2, &ss1, &ss2, 
			&nn1, &nn2, 
			&ishift, &jshift, 
			&dir, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_write_1d_hdf(char* filename, char* dsetname, 
		int& nn, int& ssp, int& nnp, int& ishift, 
		int& dir, double* phi){
	__mod_inout_MOD_write_1d_hdf(filename, dsetname,
			&nn, &ssp, &nnp, &ishift,
			&dir, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_read2_1d_hdf(char* filename, char* dsetname, 
		int& nn, int& ssp, int& nnp, int& ishift, 
		int& dir, double* phi){
	__mod_inout_MOD_read2_1d_hdf(filename, dsetname,
			&nn, &ssp, &nnp, &ishift, 
			&dir, phi, strlen(filename), strlen(dsetname));
}
inline void mpt_write_stats_hdf_4d(char* filename, char* dsetname,  
		int& ss1, int& ss2, int& ss3, int& ss4, 
		int& nn1, int& nn2, int& nn3, int& nn4, double* phi){
	__mod_inout_MOD_write_stats_hdf_4d(filename, dsetname, 
			&ss1, &ss2, &ss3, &ss4,
			&nn1, &nn2, &nn3, &nn4, 
			phi, strlen(filename), strlen(dsetname));
}
inline void mpt_read_stats_hdf_4d(char* filename, char* dsetname,
		int& ss1, int& ss2, int& ss3, int& ss4, 
		int& nn1, int& nn2, int& nn3, int& nn4, double* phi){
	__mod_inout_MOD_read_stats_hdf_4d(filename, dsetname, 
			&ss1, &ss2, &ss3, &ss4,
			&nn1, &nn2, &nn3, &nn4, 
			phi, strlen(filename), strlen(dsetname));
}
inline void mpt_start_hdf5_for_testing(){__mod_inout_MOD_start_hdf5_for_testing();}
inline void mpt_stop_hdf5_for_testing(){__mod_inout_MOD_stop_hdf5_for_testing();}
inline void mpt_filespace_props(int& vel_dir, int* stride, 
		int& s1w, int& s2w, int& s3w, 
		int& m1w, int& m2w, int& m3w,
		int& dim1, int& dim2, int& dim3){
	__mod_inout_MOD_filespace_props(&vel_dir, stride,
			&s1w, &s2w, &s3w,
			&m1w, &m2w, &m3w,
			&dim1, &dim2, &dim3);
}
inline void mpt_lambda(){__mod_inout_MOD_lambda();}
inline void mpt_lambda_2d(){__mod_inout_MOD_lambda_2d();}
inline void mpt_write_xdmf_xml(){__mod_inout_MOD_write_xdmf_xml();}
inline void mpt_write_xdmf_timecollection(){__mod_inout_MOD_write_xdmf_timecollection();}
#endif //MOD_INOUT_WRAPPER_H
