#ifndef _T1P_TILINGS_H_
#define _T1P_TILINGS_H_


typedef struct ap_abstract0_array_t {
	int n;
	ap_abstract0_t** abs;
} ap_abstract0_array_t;

ap_abstract0_array_t t1p_tilings(ap_manager_t* man, ap_abstract0_t* env1);


#endif
