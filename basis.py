from arith_ops import *
from ec import *



def find_nA_x_coord(curve, start):
    """
    Functional equivalent of C find_nA_x_coord.

    Returns:
        (xP, hint)
        where xP is an Fp2 element and hint is an int (0..127 or 0 on failure)
    """

    # assert(!fp2_is_square(&curve->A));
    assert not fp2_is_square(curve.A)

    n = start

    # initial x = n * A
    if n == 1:
        x = fp2_copy(curve.A)
    else:
        x = fp2_mul_small(curve.A, n)

    # search until x lies on the curve
    while not is_on_curve(x, curve):
        x = fp2_add(x, curve.A)
        n += 1

    # encode hint
    hint = n if n < 128 else 0
    return x, hint


def find_nqr_factor(curve, start):
    """
    Functional equivalent of C find_nqr_factor.

    Returns:
        (xP, hint)
        where xP is an Fp2 element and hint is an int (0..127 or 0 on failure)
    """

    found = False
    n = start

    while not found:
        # find b such that 1 + b^2 is a non-quadratic residue in Fp
        while True:
            tmp = n * n + 1
            if not fp_is_square(tmp):
                break
            n += 1

        b = n - 1

        # z = 1 + i*b
        z = Fp2(re=1, im=b)

        # t0 = z - 1 = i*b
        t0 = Fp2(re=0, im=b)

        # Compute A^2*(z-1) - z^2
        t1 = fp2_sqr(curve.A)
        t0 = fp2_mul(t0, t1)      # A^2 * (z - 1)
        t1 = fp2_sqr(z)
        t0 = fp2_sub(t0, t1)      # A^2 * (z - 1) - z^2

        found = not fp2_is_square(t0)

        n += 1  # mirrors the C loop behavior

    # x = -A / (1 + i*b)
    x = fp2_inv(z)
    x = fp2_mul(x, curve.A)
    x = fp2_neg(x)

    # encode hint
    hint = (n - 1) if n <= 128 else 0

    return x, hint




def ec_curve_to_basis_2f_from_hint(PQ2: ECBasis, curve: ECCurve, f: int, hint: int) -> list[bool, ECCurve, ECBasis]:
    # Normalise (A/C : 1) and ((A + 2)/4 : 1)
    
    curve = ec_normalize_curve_and_A24(curve)

    
    if fp2_is_zero(curve.A):
        PQ2 = ec_basis_E0_2f(curve, f)
        return True, curve, PQ2
    
    # Decode hint
    hint_A = hint & 1
    hint_P = hint >> 1

    # Local points
    P = PQ2.P
    Q = PQ2.Q

    if not hint_P:
        if not hint_A:
            P.x, hint = find_nA_x_coord(curve, 128)
        else:
            P.x, hint = find_nqr_factor(curve, 128)
    else:
        if not hint_A:
            P.x = fp2_mul_small(curve.A, hint_P)
        else:
            P.x.re = fp_set_one()
            P.x.im = fp_set_small(hint_P)
            P.x = fp2_inv(P.x)
            P.x = fp2_mul(P.x, curve.A)
            P.x = fp2_neg(P.x)
    

    P.z = fp2_set_one()

    # set xQ = -xP - A
    Q.x = fp2_add(curve.A, P.x)
    Q.x = fp2_neg(Q.x)
    Q.z = fp2_set_one()

    # clear cofactor
    # P = clear_cofactor_for_maximal_even_order(P, curve, f)
    # Q = clear_cofactor_for_maximal_even_order(Q, curve, f)

    P = xMUL(P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, curve)
    Q = xMUL(Q, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, curve)

    # compute PmQ
    PQ2.Q = difference_point(P, Q, curve)
    PQ2.P = copy_ec_point(P)
    PQ2.PmQ = copy_ec_point(Q)
    
    return True, curve, PQ2


def is_on_curve(x: Fp2, curve: ECCurve) -> bool:
    """
    Functional equivalent of C is_on_curve.
    Checks whether x is an x-coordinate of a point on the curve.
    """

    # C code assumes curve.C == 1
    assert fp2_is_one(curve.C)

    # t0 = x + A
    t0 = fp2_add(x, curve.A)

    # t0 = x^2 + A*x
    t0 = fp2_mul(t0, x)

    # t0 = x^2 + A*x + 1
    t0 = fp2_add_one(t0)

    # t0 = x^3 + A*x^2 + x
    t0 = fp2_mul(t0, x)

    return fp2_is_square(t0)

