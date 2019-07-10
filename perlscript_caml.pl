while (<>) {
    s/Manager.ap_manager_ptr/t Apron.Manager.t/g; 
    s/external t1p_/external /g;  
    s/Abstract0.ap_abstract0_ptr/t Apron.Abstract0.t/g;
    s/type ap_abstract0_array_t = /type t type ap_abstract0_array_t = /g;
   print;
}
