#ifndef MOD_EXCHANGE_WRAPPER_H
#define MOD_EXCHANGE_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
	void __mod_exchange_MOD_exchange(int* dir, int* vel_dir, double* phi);
	void __mod_exchange_MOD_exchange2(int* dir, int* vel_dir, int* ss1,
			int* ss2, int* ss3, int* nn1, int* nn2, int* nn3, double* phi);
	void __mod_exchange_MOD_exchange_all_all(int* vel_yes, double* phi);
	void __mod_exchange_MOD_exchange_relax(int* g, int* fb1, int* fb2, int* fb3,
			int* vel_dir, int* mirror_yes, double* phi);
	void __mod_exchange_MOD_exchange_part(int* dir, int* vel_dir, int* ss,
			double* phi);
#ifdef __cplusplus
}
#endif

// function wrappers
inline void mpt_exchange(int& dir, int& vel_dir, double* phi){
	__mod_exchange_MOD_exchange(&dir,&vel_dir, phi);
}
inline void mpt_exchange2(int& dir, int& vel_dir, int& ss1, int& ss2, int& ss3,
		int& nn1, int& nn2, int& nn3, double* phi){
	__mod_exchange_MOD_exchange2(&dir,&vel_dir, &ss1, &ss2, &ss3, &nn1, &nn2, &nn3, phi);
}
inline void mpt_exchange_all_all(int vel_yes, double* phi){
	__mod_exchange_MOD_exchange_all_all(&vel_yes, phi);
}
inline void mpt_exchange_relax(int& g, int& fb1, int& fb2, int& fb3, int& vel_dir,
		int& mirror_yes, double* phi){
	__mod_exchange_MOD_exchange_relax(&g, &fb1,&fb2, &fb3, &vel_dir, &mirror_yes, phi);
}
inline void mpt_exchange_part(int& dir, int& vel_dir, int& ss, double* phi){
	__mod_exchange_MOD_exchange_part(&dir, &vel_dir, &ss, phi);
}
#endif //MOD_EXCHANGE_WRAPPER_H
