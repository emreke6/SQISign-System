from arith_ops import *
from typing import Optional
from ec import *
from basis import *
from xisog import *
from theta_isogenies import *
from hash import *

import json

def check_canonical_basis_change_matrix(sig: Signature) -> int:
    # aux = 1 << (SQIsign_response_length + HD_extra_torsion - sig->backtracking)
    shift = SQIsign_response_length + HD_extra_torsion - sig.backtracking
    aux = 1 << shift

    for i in range(2):
        for j in range(2):
            entry = sig.mat_Bchall_can_to_B_chall[i][j]
            if mp_compare_int(aux, entry) <= 0:
                return 0

    return 1





def compute_challenge_verify(
    E_chall: ECCurve,
    sig: Signature,
    Epk: ECCurve,
    hint_pk: int
) -> bool:

    # ec_basis_t bas_EA;
    bas_EA = ECBasis(
        P=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        Q=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        PmQ=ECPoint(Fp2(0, 0), Fp2(0, 0)),
    )

    # ec_isog_even_t phi_chall;
    phi_chall = ECIsogEven(
        curve=ECCurve(
            A=Fp2(0, 0),
            C=Fp2(0, 0),
            A24=None,
            is_A24_computed_and_normalized=False,
        ),
        kernel=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        length=0,
    )

    # copy_curve(&phi_chall.curve, Epk);
    phi_chall.curve = copy_curve(Epk)

    # phi_chall.length = TORSION_EVEN_POWER - sig->backtracking;
    phi_chall.length = TORSION_EVEN_POWER - sig.backtracking

    # Compute the basis from the supplied hint
    ret1, phi_chall.curve, bas_EA = ec_curve_to_basis_2f_from_hint(
        bas_EA,
        phi_chall.curve,
        TORSION_EVEN_POWER,
        hint_pk,
    )
    if not ret1:
        return False

    # recovering the exact challenge
    ret, phi_chall.kernel = ec_ladder3pt(
        sig.chall_coeff,
        bas_EA.P,
        bas_EA.Q,
        bas_EA.PmQ,
        phi_chall.curve
    )
    if not ret:
        return False

    # Double the kernel until it has the correct order
    phi_chall.kernel,  phi_chall.curve = ec_dbl_iter(phi_chall.kernel, sig.backtracking, phi_chall.kernel, phi_chall.curve)

    # Compute the codomain
    E_chall = copy_curve(phi_chall.curve)

    ret_ec_eval, new_E_chall = ec_eval_even(E_chall, phi_chall, None, 0)

    E_chall = new_E_chall

    if ret_ec_eval:
        return False

    return True, E_chall


def ec_eval_even_strategy(
    curve_in: ECCurve,
    points: list[ECPoint],
    len_points: int,
    kernel: ECPoint,
    isog_len: int
) -> int:
    """
    In-place equivalent of the C ec_eval_even_strategy.
    Mutates `curve` and `points`.
    Returns 0 on success, -1 on failure.
    """

    # Normalize A24 in place
    curve_norm = ec_curve_normalize_A24(curve_in)
    curve = ec_curve_init()
    curve.A = curve_norm.A
    curve.C = curve_norm.C
    curve.A24 = curve_norm.A24
    curve.is_A24_computed_and_normalized = curve_norm.is_A24_computed_and_normalized

    # Local copy of A24
    A24 = ec_point_init()
    A24 = copy_ec_point(curve.A24)

    # Compute stack space
    space = 1
    i = 1
    while i < isog_len:
        space += 1
        i *= 2

    # Stack of kernel splits and remaining orders
    splits = [ec_point_init() for i in range(space)]
    todo = [0 for i in range(space)] 

    splits[0] = copy_ec_point(kernel)
    todo[0] = isog_len
    current = 0

    # Chain of 4-isogenies
    for j in range(isog_len // 2):
        assert current >= 0
        assert todo[current] >= 1

        # Find next point of order 4
        ctr = 0
        while todo[current] != 2:
            assert todo[current] >= 3
            current += 1
            assert current < space

            splits[current] = copy_ec_point(splits[current - 1])

            num_dbls = (todo[current - 1] // 4) * 2 + (todo[current - 1] % 2)
            todo[current] = todo[current - 1] - num_dbls

            
            while num_dbls > 0:
                splits[current] = xDBL_A24(splits[current], A24, False)
                num_dbls -= 1

        if j == 0:
            assert fp2_is_one(A24.z)
            if not ec_is_four_torsion(splits[current], curve):
                return -1
            
            T = xDBL_A24(splits[current], A24, False)
            if fp2_is_zero(T.x):
                return -1
        else:
            assert todo[current] == 2

        kps4 = ECKPS4()
        # Evaluate 4-isogeny
        kps4, A24  = xisog_4(A24, splits[current])

        # Update splits in place
        new_splits = xeval_4(splits, current, kps4)
        for i in range(current):
            splits[i] = new_splits[i]

        for i in range(current):
            todo[i] -= 2

        # Update points in place
        # points = xeval_4(points, len_points, kps4)
        new_points = xeval_4(points, len_points, kps4)
        for i in range(len_points):
            points[i] = new_points[i]

        current -= 1

    # Final 2-isogeny if needed
    if isog_len % 2:
        print("isog_len % 2")
        if isog_len == 1 and not ec_is_two_torsion(splits[0], curve):
            return -1

        if fp2_is_zero(splits[0].x):
            return -1

        kps2, A24 = xisog_2(A24, splits[0])
        updated_points = xeval_2(points, len_points, kps2)
        for i in range(len_points):
            points[i] = updated_points[i]

    # Convert A24 back to (A:C) in place
    curve_ac = A24_to_AC(curve, A24)
    curve.A = curve_ac.A
    curve.C = curve_ac.C
    curve.A24 = curve_ac.A24
    curve.is_A24_computed_and_normalized = False

    return 0, curve

def ec_eval_even(
    image: ECCurve,
    phi: ECIsogEven,
    points: list[ECPoint],
    len_points: int
) -> int:
    """
    In-place equivalent of C ec_eval_even.
    """

    curve_copy = copy_curve(phi.curve)
    image.A = curve_copy.A
    image.C = curve_copy.C
    image.A24 = curve_copy.A24
    image.is_A24_computed_and_normalized = curve_copy.is_A24_computed_and_normalized

    ret1, new_curve = ec_eval_even_strategy(
        image,
        points,
        len_points,
        phi.kernel,
        phi.length
    )

    image = new_curve

    return ret1, new_curve

def mp_sub(a: int, b: int):
    return a-b

def mp_mod_2exp(a: int, e: int) -> int:
    return a & ((1 << e) - 1)

def matrix_scalar_application_even_basis(bas: ECBasis,
                                         E: ECCurve,
                                         mat: list,
                                         f: int) -> int:
    # scalar_t scalar0, scalar1
    scalar0 = 0
    scalar1 = 0

    tmp_bas = ECBasis(
        P=ec_point_init(),
        Q=ec_point_init(),
        PmQ=ec_point_init()
    )
    tmp_bas.P = bas.P
    tmp_bas.Q = bas.Q
    tmp_bas.PmQ = bas.PmQ

    # R = [a]P + [b]Q
    bas.P, ret = ec_biscalar_mul(
        mat[0][0],
        mat[1][0],
        f,
        tmp_bas,
        E
    )
    if not ret:
        return 0

    # S = [c]P + [d]Q
    bas.Q, ret2 = ec_biscalar_mul(
        mat[0][1],
        mat[1][1],
        f,
        tmp_bas,
        E
    )
    if not ret2:
        return 0

    # scalar0 = a - c mod 2^f
    scalar0 = mp_sub(mat[0][0], mat[0][1])
    scalar0 = mp_mod_2exp(scalar0, f)

    # scalar1 = b - d mod 2^f
    scalar1 = mp_sub(mat[1][0], mat[1][1])
    scalar1 = mp_mod_2exp(scalar1, f)


    bas.PmQ, ret3 = ec_biscalar_mul(
        scalar0,
        scalar1,
        f,
        tmp_bas,
        E
    )
    if not ret3:
        return 0

    return bas, 1



def challenge_and_aux_basis_verify(B_chall_can: ECBasis, B_aux_can: ECBasis, E_chall: ECCurve, E_aux: ECCurve, sig: Signature, pow_dim2_deg_resp: int):
    


    ret1, E_chall ,B_chall_can = ec_curve_to_basis_2f_from_hint(B_chall_can, E_chall, TORSION_EVEN_POWER, sig.hint_chall)
 
    if not ret1:
        return False
    
    B_chall_can, E_chall = ec_dbl_iter_basis(B_chall_can, 
                                             TORSION_EVEN_POWER - pow_dim2_deg_resp - HD_extra_torsion - sig.two_resp_length,
                                             B_chall_can,
                                             E_chall)

    ret1, E_aux, B_aux_can = ec_curve_to_basis_2f_from_hint(B_aux_can, E_aux, TORSION_EVEN_POWER, sig.hint_aux)

    B_aux_can, E_aux = ec_dbl_iter_basis(B_aux_can, TORSION_EVEN_POWER - pow_dim2_deg_resp - HD_extra_torsion, B_aux_can, E_aux)


    B_chall_can, ret1 = matrix_scalar_application_even_basis(B_chall_can,
                                                E_chall,
                                                sig.mat_Bchall_can_to_B_chall,
                                                pow_dim2_deg_resp + HD_extra_torsion + sig.two_resp_length)
    
    return B_chall_can, B_aux_can , E_aux, E_chall ,True


def compute_commitment_curve_verify(E_com: ECCurve,
                                B_chall_can: ECBasis,
                                B_aux_can: ECBasis,
                                E_chall: ECCurve,
                                E_aux: ECCurve,
                                pow_dim2_deg_resp: int):
    

    EchallxEaux = theta_couple_curve_t(E1=ec_curve_init(), E2=ec_curve_init())
    EchallxEaux.E1 = E_chall
    EchallxEaux.E2 = E_aux

    dim_two_ker = copy_bases_to_kernel(B_chall_can, B_aux_can)

    codomain = theta_couple_curve_t(E1=ec_curve_init(), E2=ec_curve_init())
    codomain.E1 = ec_curve_init()
    codomain.E2 = ec_curve_init()


    if pow_dim2_deg_resp == 0:
        codomain_splits = 1
        codomain.E1 = copy_curve(EchallxEaux.E1)
        codomain.E2 = copy_curve(EchallxEaux.E2)
        # We still need to check that E_chall is supersingular
        # This assumes that HD_extra_torsion == 2
        if not ec_is_basis_four_torsion(B_chall_can, E_chall):
            return 0
    else:
        codomain_splits, EchallxEaux, codomain  = theta_chain_compute_and_eval_verify(
            pow_dim2_deg_resp, EchallxEaux, dim_two_ker, True, 0)
    

    # computing the commitment curve
    # its always the first one because of our (2^n,2^n)-isogeny formulae
    E_com = copy_curve(codomain.E1)

    return codomain_splits, E_com


def mp_is_odd(x: int) -> bool:
    return (x & 1) == 1

def mp_is_even(x: int) -> bool:
    return (x & 1) == 0


# naive implementation
def ec_eval_small_chain(curve: ECCurve,
                        kernel: ECPoint,
                        len: int,
                        points_in: list[ECPoint],
                        len_points: int,
                        special: bool): #  do we allow special isogenies?
    curve_new = ECCurve(A=curve.A, C = curve.C, A24=curve.A24, is_A24_computed_and_normalized=curve.is_A24_computed_and_normalized)
    # curve_new.A = curve.A
    # curve_new.C = curve.C
    # curve_new.A24 = curve.A24
    # curve_new.is_A24_computed_and_normalized = curve.is_A24_computed_and_normalized

    points = [points_in[i] for i in range(len_points)]
    A24 = ec_point_init()
    A24 = AC_to_A24(curve_new)

    kps= ECKPS2()
    small_K = ec_point_init()
    big_K = ec_point_init()
    big_K = copy_ec_point(kernel)

    for i in range(len):
        small_K = copy_ec_point(big_K)
        # small_K = big_K;
        for j in range(len - i - 1):
            print("dbl", i, j)
            small_K = xDBL_A24(small_K, A24, False)

        # Check the order of the point before the first isogeny step
        if i == 0 and not ec_is_two_torsion(small_K, curve_new):
            return -1
        # Perform isogeny step
        if fp2_is_zero(small_K.x):
            if special:
                B24 = ec_point_init()
                kps = xisog_2_singular(kps, B24, A24)
                big_K = xeval_2_singular(big_K, 1, kps)
                points = xeval_2_singular(points, len_points, kps)
                A24 = copy_ec_point(B24)
            else:
                return -1
            
        else:
            kps, A24    = xisog_2(A24, small_K)
            big_K       = xeval_2_one(big_K, 1, kps)
            points      = xeval_2(points, len_points, kps)
        

 
    curve_new = A24_to_AC(curve_new, A24)

    curve_new.is_A24_computed_and_normalized = False
    return 0, curve_new, points




def two_response_isogeny_verify(E_chall: ECCurve,
                                B_chall_can: ECBasis,
                                sig: Signature,
                                pow_dim2_deg_resp: int) -> int:
    """
    Verifies the 2-isogeny response and updates B_chall_can in place.
    Returns 1 on success, 0 on failure.
    """

    ker = ec_point_init()
    points = [ec_point_init() for _ in range(3)]

    # Choose kernel for small 2-isogenies
    if (mp_is_even(sig.mat_Bchall_can_to_B_chall[0][0]) and
        mp_is_even(sig.mat_Bchall_can_to_B_chall[1][0])):
        ker = copy_ec_point(B_chall_can.Q)
    else:
        ker = copy_ec_point(B_chall_can.P)

    # Copy basis points
    points[0] = copy_ec_point(B_chall_can.P)
    points[1] = copy_ec_point(B_chall_can.Q)
    points[2] = copy_ec_point(B_chall_can.PmQ)

    # Multiply kernel by 2^(pow_dim2_deg_resp + HD_extra_torsion)
    ker, E_chall = ec_dbl_iter(
        ker,
        pow_dim2_deg_resp + HD_extra_torsion,
        ker,
        E_chall
    )


    # Evaluate the small 2-isogeny chain
    ret, E_chall, points,  = ec_eval_small_chain(
        E_chall,
        ker,
        sig.two_resp_length,
        points,
        3,
        False
    )
    if ret:
        return 0
    

    # Update canonical challenge basis
    B_chall_can.P   = copy_ec_point(points[0])
    B_chall_can.Q   = copy_ec_point(points[1])
    B_chall_can.PmQ = copy_ec_point(points[2])

    return E_chall, B_chall_can, 1



def mp_compare(a: int, b: int):
    # Multiprecision comparison, a=b? : (1) a>b, (0) a=b, (-1) a<b

    
    if (a > b):
        return 1
    elif a < b:
        return -1

    return 0


def protocols_verify(sig: Signature, pk: PublicKey, m: bytes, l: int, orig_curve: ECCurve) -> bool:

    verify = 0

    if not check_canonical_basis_change_matrix(sig):
        return False

    pow_dim2_deg_resp = SQIsign_response_length - sig.two_resp_length - sig.backtracking

    if pow_dim2_deg_resp < 0 or pow_dim2_deg_resp == 1:
        return False

    if not ec_curve_verify_A(pk.curve.A):
        return False

    E_aux = ec_curve_init_from_A(sig.E_aux_A)

    if E_aux is None:
        return False
    
    assert fp2_is_one(pk.curve.C)
    assert not pk.curve.is_A24_computed_and_normalized

    E_chall = ec_curve_init()
    ret_cmpt, E_chall_new = compute_challenge_verify(E_chall, sig_obj, pk.curve, pk.hint_pk)


    E_chall = E_chall_new

    if not ret_cmpt:
        return False
    

    B_aux_can = ECBasis(
        P=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        Q=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        PmQ=ECPoint(Fp2(0, 0), Fp2(0, 0)),
    )
    B_chall_can =  ECBasis(
        P=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        Q=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        PmQ=ECPoint(Fp2(0, 0), Fp2(0, 0)),
    )

    B_chall_can, B_aux_can , E_aux, E_chall, ret_chal_verify = challenge_and_aux_basis_verify(
            B_chall_can, B_aux_can,
            E_chall, E_aux,
            sig, pow_dim2_deg_resp)
    
    if not ret_chal_verify:
        return False


    if sig.two_resp_length > 0:
        E_chall, B_chall_can, ret1 = two_response_isogeny_verify(E_chall, B_chall_can,
                                                                    sig, pow_dim2_deg_resp)
        
        if not ret1:
            return 0


    E_com = ec_curve_init()
    ret_cmt_vrf, E_com = compute_commitment_curve_verify(
            E_com,
            B_chall_can, B_aux_can,
            E_chall, E_aux,
            pow_dim2_deg_resp)
        
    if not ret_cmt_vrf:    
        return 0

    m_bytes = m.to_bytes(length, "little")
    orig_chall = hash_to_challenge(pk_obj, orig_curve  , m_bytes, l,  FP2_ENCODED_BYTES, SECURITY_BITS, HASH_ITERATIONS, TORSION_EVEN_POWER, SQIsign_response_length, RADIX)
    
    chk_chall = hash_to_challenge(pk_obj, E_com, m_bytes, l, FP2_ENCODED_BYTES, SECURITY_BITS, HASH_ITERATIONS, TORSION_EVEN_POWER, SQIsign_response_length, RADIX)

    verify = mp_compare(orig_chall, chk_chall) == 0
    return verify



if __name__ == "__main__":
    # ---- load parsed JSON ----
    with open("parsed.json") as f:
        data = json.load(f)

    signatures_arr = data["signatures"]
    public_key_arr = data["public_keys"]
    Ecom_Eaux_E1_arr = data["Ecom_Eaux_E1"]

    # convert messages to your previous format
    message_arr = [int(m["message"], 16) for m in data["messages"]]
    length_arr = [m["length"] for m in data["messages"]]
    # ---- your verification loop ----
    # test_size = len(message_arr)

    test_size = 1

    for i in range(test_size):
        sig_obj = Signature(
            E_aux_A=Fp2(
                re=signatures_arr[i]["E_aux_A"]["real"],
                im=signatures_arr[i]["E_aux_A"]["imag"]
            ),
            mat_Bchall_can_to_B_chall=signatures_arr[i]["mat_Bchall_can_to_B_chall"],
            backtracking=signatures_arr[i]["backtracking"],
            two_resp_length=signatures_arr[i]["two_resp_length"],
            hint_aux=signatures_arr[i]["hint_aux"],
            hint_chall=signatures_arr[i]["hint_chall"],
            chall_coeff=signatures_arr[i]["chall_coeff"]
        )

        pk_obj = PublicKey(
            curve=ECCurve(
                A=Fp2(
                    public_key_arr[i]["curve"]["A"]["real"],
                    public_key_arr[i]["curve"]["A"]["imag"]
                ),
                C=Fp2(
                    public_key_arr[i]["curve"]["C"]["real"],
                    public_key_arr[i]["curve"]["C"]["imag"]
                ),
                A24=None,
                is_A24_computed_and_normalized=False
            ),
            hint_pk=public_key_arr[i]["hint_pk"]
        )

        orig_curve_obj = ECCurve(
            A=Fp2(
                Ecom_Eaux_E1_arr[i]["curve"]["A"]["real"],
                Ecom_Eaux_E1_arr[i]["curve"]["A"]["imag"]
            ),
            C=Fp2(
                Ecom_Eaux_E1_arr[i]["curve"]["C"]["real"],
                Ecom_Eaux_E1_arr[i]["curve"]["C"]["imag"]
            ),
            A24=ECPoint(
                x=Fp2(
                    Ecom_Eaux_E1_arr[i]["curve"]["A24"]["x"]["real"],
                    Ecom_Eaux_E1_arr[i]["curve"]["A24"]["x"]["imag"]
                ),
                z=Fp2(
                    Ecom_Eaux_E1_arr[i]["curve"]["A24"]["z"]["real"],
                    Ecom_Eaux_E1_arr[i]["curve"]["A24"]["z"]["imag"]
                ),
            ),
            is_A24_computed_and_normalized=
                Ecom_Eaux_E1_arr[i]["curve"]["is_A24_computed_and_normalized"]
        )

        message = message_arr[i]
        length = length_arr[i]
                

        result = protocols_verify(sig_obj, pk_obj, message, length, orig_curve_obj)
        
        print("protocols_verify =", i ,result)


    