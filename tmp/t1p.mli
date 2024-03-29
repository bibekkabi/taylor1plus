(* File generated from t1p.idl *)

type ap_abstract0_array_t = Abstract0.ap_abstract0_ptr array

(** Taylor1+ abstract domain (beta version) *)

(** Type of Taylor1+ forms.

Each dimension/variable [x_i] has the affine form: alpha_0 + Sum(alpha_i*eps_i), where eps_i are the noise symbols, and alpha_i their associated coefficients.

Abstract values which are Taylor1+ forms (affine forms) have the type [t Apron.AbstractX.t].

Managers allocated for Taylor1+ abstract values have the type [t Apron.manager.t].
*)


(** Create a Taylor1+ manager. *)
external t1p_manager_alloc : unit -> Manager.ap_manager_ptr
	= "camlidl_t1p_t1p_manager_alloc"



(**
{2 Compilation information}

{3 Bytecode compilation}
To compile to bytecode, you should first generate a custom
interpreter with a command which should look like:

[ocamlc -I $APRON_PREFIX/lib -make-runtime -o myrun bigarray.cma gmp.cma apron.cma t1p.cma]

and then you compile and link your example [X.ml] with

[ocamlc -I $APRON_PREFIX/lib -c X.ml] and

[ocamlc -I $APRON_PREFIX/lib -use-runtime myrun -o X bigarray.cma gmp.cma apron.cma t1p.cma X.cmo]

{b Comments:} The C libraries related to [gmp.cma] and [apron.cma] are
automatically looked for (thanks to the auto-linking feature provided by
[ocamlc]). For [t1p.cma], the library [libt1p.a], identic to [libt1pMPQ.a], is
selected by default. The [-noautolink] option should be used to select a differetn version. See the C documentation of [t1p] library for details.

With the [-noautolink] option, the generation of the custom
runtime executable should be done with

[ocamlc -I $APRON_PREFIX/lib -noautolink -make-runtime -o myrun bigarray.cma gmp.cma apron.cma t1p.cma -ccopt "-L$GMP_PREFIX/lib ..." -cclib "-lt1p_caml -lt1pMPQ -lapron_caml -lapron -lgmp_caml -lmpfr -lgmp -lbigarray -lcamlidl"]

{3 Native-code compilation}
You compile and link with

[ocamlopt -I $APRON_PREFIX/lib -c X.ml] and

[ocamlopt -I $APRON_PREFIX/lib -o X bigarray.cmxa gmp.cmxa apron.cmxa t1p.cmxa X.cmx]

{b Comments:} Same as for bytecode compilation. With the
[-noautolink] option, the linking command becomes

[ocamlopt -I $APRON_PREFIX/lib -o X bigarray.cmxa gmp.cmxa apron.cmxa t1p.cmxa -ccopt "-L$GMP_PREFIX/lib ..." -cclib "-lt1p_caml -lt1pMPQ -lapron_caml -lapron -lgmp_caml -lmpfr -lgmp -lbigarray -lcamlidl" X.cmx]
*)

(** Tiling. *)
external t1p_tilings : Manager.ap_manager_ptr -> Abstract0.ap_abstract0_ptr -> ap_abstract0_array_t
	= "camlidl_t1p_t1p_tilings"

