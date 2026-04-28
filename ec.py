from arith_ops import *
from typing import Optional
from ec import *



def ec_point_init() -> ECPoint:
    return ECPoint(
        x=fp2_set_one(),
        z=fp2_set_zero(),
    )

def ec_curve_init() -> ECCurve:
    return ECCurve(
        A=fp2_set_zero(),
        C=fp2_set_one(),
        A24=ec_point_init(),
        is_A24_computed_and_normalized=False,
    )


@dataclass
class ECBasis:
    P: ECPoint
    Q: ECPoint
    PmQ: ECPoint


@dataclass
class ECIsogEven:
    curve: ECCurve
    kernel: ECPoint
    length: int

@dataclass
class ECKPS2:
    K: ECPoint = field(
        default_factory=lambda: ECPoint(Fp2(0, 0), Fp2(0, 0))
    )

@dataclass
class ECKPS4:
    K: list[ECPoint] = field(default_factory=lambda: [
        ECPoint(Fp2(0, 0), Fp2(0, 0)),
        ECPoint(Fp2(0, 0), Fp2(0, 0)),
        ECPoint(Fp2(0, 0), Fp2(0, 0)),
    ])

@dataclass
class add_components_t:
    u: Fp2 = None
    v: Fp2 = None
    w: Fp2 = None

@dataclass
class theta_couple_curve_t:
    E1: ECCurve
    E2: ECCurve

@dataclass
class theta_couple_point_t:
    P1: ECPoint = None
    P2: ECPoint = None

@dataclass
class theta_kernel_couple_points_t:
    T1: theta_couple_point_t
    T2: theta_couple_point_t
    T1m2: theta_couple_point_t

@dataclass
class jac_point_t:
    x: Fp2 = None
    y: Fp2 = None
    z: Fp2 = None

    def __post_init__(self):
        if self.x is None:
            self.x = fp2_set_zero()
        if self.y is None:
            self.y = fp2_set_zero()
        if self.z is None:
            self.z = fp2_set_zero()


def copy_jac_point(jac_point_in : jac_point_t):
    jac_new = jac_point_t(x= jac_point_in.x, y=jac_point_in.y, z=jac_point_in.z)

    return jac_new


@dataclass
class theta_couple_jac_point_t:
    P1: jac_point_t = None
    P2: jac_point_t = None

@dataclass
class theta_point_compact_t:
    x: Fp2
    y: Fp2

@dataclass
class basis_change_matrix_t:
    m: list[list[Fp2]] = field(
        default_factory=lambda: [[Fp2(0, 0) for _ in range(4)] for _ in range(4)]
    )

def basis_change_matrix_init():
    return basis_change_matrix_t(
        m=list(
            list(Fp2() for _ in range(4))
            for _ in range(4)
        )
    )

def basis_change_matrix_init_zero():
    return basis_change_matrix_t(
        m=list(
            list(Fp2(re=0, im=0) for _ in range(4))
            for _ in range(4)
        )
    )
@dataclass
class translation_matrix_t:
    g00: Fp2 = None
    g01: Fp2 = None
    g10: Fp2 = None
    g11: Fp2 = None


@dataclass
class theta_gluing_t:
    domain = theta_couple_curve_t
    xyK1_8 = theta_couple_jac_point_t
    imageK1_8 = theta_point_compact_t
    M = basis_change_matrix_t
    precomputation = ThetaPoint
    codomain = ThetaPoint

@dataclass
class theta_isogeny_t:
    T1_8: ThetaPoint = None
    T2_8: ThetaPoint = None
    hadamard_bool_1: bool = None
    hadamard_bool_2: bool = None
    domain: ThetaStructure = None
    precomputation: ThetaPoint = None
    codomain: ThetaStructure = None

@dataclass
class theta_splitting_t:
    M: basis_change_matrix_t = field(default_factory=basis_change_matrix_t)
    B: ThetaStructure = field(default_factory=ThetaStructure)


EVEN_INDEX = [
    (0, 0), (0, 1), (0, 2), (0, 3),
    (1, 0), (1, 2),
    (2, 0), (2, 1),
    (3, 0), (3, 3),
]

CHI_EVAL = [
    [ 1,  1,  1,  1],
    [ 1, -1,  1, -1],
    [ 1,  1, -1, -1],
    [ 1, -1, -1,  1],
]

FP2_ZERO = 0
FP2_ONE =  1
FP2_I = 2
FP2_MINUS_ONE = 3
FP2_MINUS_I = 4



SPLITTING_TRANSFORMS = (
    (
        (FP2_ONE,  FP2_I,        FP2_ONE,        FP2_I),
        (FP2_ONE,  FP2_MINUS_I,  FP2_MINUS_ONE,  FP2_I),
        (FP2_ONE,  FP2_I,        FP2_MINUS_ONE,  FP2_MINUS_I),
        (FP2_MINUS_ONE, FP2_I,   FP2_MINUS_ONE,  FP2_I),
    ),
    (
        (FP2_ONE,  FP2_ZERO, FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE),
        (FP2_ZERO, FP2_ZERO, FP2_ONE,  FP2_ZERO),
        (FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO, FP2_ZERO),
    ),
    (
        (FP2_ONE,  FP2_ZERO, FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ONE,  FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE),
        (FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE, FP2_ZERO),
    ),
    (
        (FP2_ONE,  FP2_ZERO, FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ONE,  FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ONE,  FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_MINUS_ONE),
    ),
    (
        (FP2_ONE,  FP2_ONE,        FP2_ONE,        FP2_ONE),
        (FP2_ONE,  FP2_MINUS_ONE,  FP2_MINUS_ONE,  FP2_ONE),
        (FP2_ONE,  FP2_ONE,        FP2_MINUS_ONE,  FP2_MINUS_ONE),
        (FP2_MINUS_ONE, FP2_ONE,   FP2_MINUS_ONE,  FP2_ONE),
    ),
    (
        (FP2_ONE,  FP2_ZERO, FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ONE,  FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE),
        (FP2_ZERO, FP2_ZERO, FP2_ONE,  FP2_ZERO),
    ),
    (
        (FP2_ONE,  FP2_ONE,        FP2_ONE,        FP2_ONE),
        (FP2_ONE,  FP2_MINUS_ONE,  FP2_ONE,        FP2_MINUS_ONE),
        (FP2_ONE,  FP2_MINUS_ONE,  FP2_MINUS_ONE,  FP2_ONE),
        (FP2_MINUS_ONE, FP2_MINUS_ONE, FP2_ONE, FP2_ONE),
    ),
    (
        (FP2_ONE,  FP2_ONE,        FP2_ONE,        FP2_ONE),
        (FP2_ONE,  FP2_MINUS_ONE,  FP2_ONE,        FP2_MINUS_ONE),
        (FP2_ONE,  FP2_MINUS_ONE,  FP2_MINUS_ONE,  FP2_ONE),
        (FP2_ONE,  FP2_ONE,        FP2_MINUS_ONE,  FP2_MINUS_ONE),
    ),
    (
        (FP2_ONE,  FP2_ONE,        FP2_ONE,        FP2_ONE),
        (FP2_ONE,  FP2_MINUS_ONE,  FP2_ONE,        FP2_MINUS_ONE),
        (FP2_ONE,  FP2_ONE,        FP2_MINUS_ONE,  FP2_MINUS_ONE),
        (FP2_MINUS_ONE, FP2_ONE,   FP2_ONE,        FP2_MINUS_ONE),
    ),
    (
        (FP2_ONE,  FP2_ZERO, FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ONE,  FP2_ZERO, FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ONE,  FP2_ZERO),
        (FP2_ZERO, FP2_ZERO, FP2_ZERO, FP2_ONE),
    ),
)

FP2_CONSTANTS = (
    # 0
    Fp2(0, 0),

    # 1
    Fp2(
        0x33
        + (0x0 << 64)
        + (0x0 << 128)
        + (0x100000000000000 << 192),
        0
    ),
    # i
    Fp2(
        0,
        0x33
        + (0x0 << 64)
        + (0x0 << 128)
        + (0x100000000000000 << 192)
    ),

    # -1
    Fp2(
        0xffffffffffffffcc
        + (0xffffffffffffffff << 64)
        + (0xffffffffffffffff << 128)
        + (0x3ffffffffffffff << 192),
        0
    ),

    # -i
    Fp2(
        0,
        0xffffffffffffffcc
        + (0xffffffffffffffff << 64)
        + (0xffffffffffffffff << 128)
        + (0x3ffffffffffffff << 192)
    ),
)


def copy_bases_to_kernel(B1: ECBasis, B2: ECBasis)->theta_kernel_couple_points_t:

    ker = theta_kernel_couple_points_t(
        T1 = theta_couple_point_t(P1=ec_point_init(), P2=ec_point_init),
        T2 = theta_couple_point_t(P1=ec_point_init(), P2=ec_point_init),
        T1m2 = theta_couple_point_t(P1=ec_point_init(), P2=ec_point_init)
    )
    #Copy the basis on E1 to (P, _) on T1, T2 and T1 - T2
    ker.T1.P1 = copy_ec_point(B1.P)
    ker.T2.P1 = copy_ec_point(B1.Q)
    ker.T1m2.P1 = copy_ec_point(B1.PmQ)

    #Copy the basis on E2 to (_, P) on T1, T2 and T1 - T2
    ker.T1.P2 = copy_ec_point(B2.P)
    ker.T2.P2 = copy_ec_point(B2.Q)
    ker.T1m2.P2 = copy_ec_point(B2.PmQ)

    return ker


def multiple_mp_shiftl_int(x: int, shift: int) -> int:
    """Equivalent of multiple_mp_shiftl but using Python int"""
    return x << shift


def mp_compare_int(a: int, b: int) -> int:
    if a > b: return 1
    if a < b: return -1
    return 0





def A24_to_AC(E: ECCurve, A24: ECPoint)->ECCurve:
    new_curve = ECCurve(A=E.A, C=E.C, A24=E.A24, is_A24_computed_and_normalized=E.is_A24_computed_and_normalized)

    aa1 = fp2_add(A24.x, A24.x)
    aa1 = fp2_sub(aa1, A24.z)
    aa1 = fp2_add(aa1, aa1)
    new_curve.A = aa1

    new_curve.C = fp2_copy(A24.z)

    return new_curve


def ec_curve_verify_A(A: Fp2) -> bool:
    # t = 2
    t = Fp2(2, 0)
    if fp2_is_equal(A, t):
        return False

    # t = -2
    t = Fp2(-2, 0)
    if fp2_is_equal(A, t):
        return False

    return True

def ec_curve_init_from_A(A: Fp2) -> ECCurve:
    E = ec_curve_init()
    E.A = fp2_copy(A)

    if not ec_curve_verify_A(A):
        return None

    return E




def ec_curve_verify_A(A: Fp2) -> bool:
    """
    Verify the Montgomery coefficient A is valid (A^2 - 4 != 0)
    Equivalent C logic:
      - reject if A ==  2
      - reject if A == -2
      - otherwise accept
    """

    # t = 1
    t = fp2_set_one()

    # t = 2
    t = fp2_add(t, t)

    # if A == 2: reject
    if fp2_is_equal(A, t):
        return False

    # t = -2
    t = Fp2(re=fp_neg(t.re), im=t.im)

    # if A == -2: reject
    if fp2_is_equal(A, t):
        return False

    return True


def copy_curve(src: ECCurve) -> ECCurve:
    return ECCurve(
        A=fp2_copy(src.A),
        C=fp2_copy(src.C),
        A24=src.A24,  # safe to share if immutable; otherwise copy_point(src.A24)
        is_A24_computed_and_normalized=src.is_A24_computed_and_normalized,
    )


def ec_normalize_curve_and_A24(E: ECCurve) -> ECCurve:
    # Step 1: normalize curve if needed
    
    new_ret = E
    if not fp2_is_one(new_ret.C):
        new_ret = ec_normalize_curve(new_ret)

    # Step 2: compute A24 if needed
    if not new_ret.is_A24_computed_and_normalized:
        # A24.x = (A + 2) / 4
        x = fp2_add_one(new_ret.A)
        x = fp2_add_one(x)

        x = fp2_half(x)
        x = fp2_half(x)

        A24 = ECPoint(
            x=x,
            z=fp2_set_one(),
        )

        E_new = ec_curve_init()

        E_new.A = new_ret.A
        E_new.C = new_ret.C
        E_new.A24 = A24
        E_new.is_A24_computed_and_normalized = True
    
        return E_new
    else:
        return new_ret

def ec_normalize_curve(E: ECCurve) -> ECCurve:
    C_inv = fp2_inv(E.C)
    A_new = fp2_mul(E.A, C_inv)

    return ECCurve(
        A=A_new,
        C=fp2_set_one(),
        A24=E.A24,
        is_A24_computed_and_normalized=E.is_A24_computed_and_normalized,
    )


def xDBL_E0(P):
    """
    Doubling of a Montgomery point on E0 with (A:C) = (0:1).
    Functional version: returns a new ECPoint.
    """

    # t0 = (XP + ZP)^2
    t0 = fp2_add(P.x, P.z)
    t0 = fp2_sqr(t0)

    # t1 = (XP - ZP)^2
    t1 = fp2_sub(P.x, P.z)
    t1 = fp2_sqr(t1)

    # t2 = t0 - t1
    t2 = fp2_sub(t0, t1)

    # t1 = 2 * t1
    t1 = fp2_add(t1, t1)

    # XQ = t0 * t1
    XQ = fp2_mul(t0, t1)

    # ZQ = (t1 + t2) * t2
    ZQ = fp2_add(t1, t2)
    ZQ = fp2_mul(ZQ, t2)

    return ECPoint(x=XQ, z=ZQ)




def copy_point(Q):
    """
    Functional equivalent of C copy_point.
    Returns a new ECPoint copied from Q.
    """
    new = Fp2(re=0, im=0)
    new.x = fp2_copy(Q.x)
    new.z = fp2_copy(Q.z)

    return new

def copy_ec_point(Q: ECPoint)-> ECPoint:
    """
    Functional equivalent of C copy_point.
    Returns a new ECPoint copied from Q.
    """
    new = ec_point_init()
    new.x = Q.x
    new.z = Q.z

    return new

def difference_point(P, Q, curve):
    """
    Functional equivalent of C difference_point.
    Returns a new ECPoint representing (P - Q).
    """



    # Temporary variables
    t0 = fp2_mul(P.x, Q.x)
    t1 = fp2_mul(P.z, Q.z)


    Bxx = fp2_sub(t0, t1)
    Bxx = fp2_sqr(Bxx)
    Bxx = fp2_mul(Bxx, curve.C)

    Bxz = fp2_add(t0, t1)

    t0 = fp2_mul(P.x, Q.z)
    t1 = fp2_mul(P.z, Q.x)

    Bzz = fp2_add(t0, t1)
    Bxz = fp2_mul(Bxz, Bzz)


    Bzz = fp2_sub(t0, t1)
    Bzz = fp2_sqr(Bzz)
    Bzz = fp2_mul(Bzz, curve.C)

    Bxz = fp2_mul(Bxz, curve.C)

    t0 = fp2_mul(t0, t1)
    t0 = fp2_mul(t0, curve.A)
    t0 = fp2_add(t0, t0)
    Bxz = fp2_add(Bxz, t0)

    

    # Normalization factor
    t0 = Fp2(curve.C.re, fp_neg(curve.C.im))
    t0 = fp2_sqr(t0)
    t0 = fp2_mul(t0, curve.C)

    t1 = Fp2(P.z.re, fp_neg(P.z.im))
    t1 = fp2_sqr(t1)
    t0 = fp2_mul(t0, t1)

    t1 = Fp2(Q.z.re, fp_neg(Q.z.im))
    t1 = fp2_sqr(t1)
    t0 = fp2_mul(t0, t1)

    Bxx = fp2_mul(Bxx, t0)
    Bxz = fp2_mul(Bxz, t0)
    Bzz = fp2_mul(Bzz, t0)

    # Solve quadratic
    t0 = fp2_sqr(Bxz)
    t1 = fp2_mul(Bxx, Bzz)
    t0 = fp2_sub(t0, t1)

    return ECPoint(x=t0, z=t1)
    
    t0 =  fp2_sqrt(t0)
    X = fp2_add(Bxz, t0)
    Z = fp2_copy(Bzz)

    return ECPoint(x=X, z=Z)



def ec_basis_E0_2f(curve:ECCurve, f: int) -> ECBasis:
    # assert(fp2_is_zero(&curve->A));
    assert fp2_is_zero(curve.A)

    # local points
    P = ECPoint(x=None, z=None)
    Q = ECPoint(x=None, z=None)

    # Set P, Q to precomputed (X : 1) values
    fp2_copy(P.x, BASIS_E0_PX)
    fp2_copy(Q.x, BASIS_E0_QX)
    fp2_set_one(P.z)
    fp2_set_one(Q.z)

    # clear the power of two to get a point of order 2^f
    for _ in range(TORSION_EVEN_POWER - f):
        P = xDBL_E0(P)
        Q = xDBL_E0(Q)

    # Set P, Q in the basis and compute x(P - Q)
    PQ2 = ECBasis(
        P=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        Q=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        PmQ=ECPoint(Fp2(0, 0), Fp2(0, 0)),
    )
    PQ2.P = copy_ec_point(P)
    PQ2.Q = copy_ec_point(Q)
    PQ2.PmQ = difference_point(P, Q, curve)

    return PQ2




def clear_cofactor_for_maximal_even_order(
    P: ECPoint,
    curve: ECCurve,
    f: int
) -> ECPoint:
    """
    Pure functional version.
    Does NOT mutate P.
    Returns a new ECPoint.
    """

    # Clear the odd cofactor
    P2 = ec_mul(
        p_cofactor_for_2f,
        P_COFACTOR_FOR_2F_BITLENGTH,
        P,
        curve
    )

    # Clear remaining powers of two
    # for _ in range(TORSION_EVEN_POWER - f):
    #     P2 = xDBL_A24(
    #         P2,
    #         curve.A24,
    #         curve.is_A24_computed_and_normalized
    #     )

    #P2 = ec_point_init()

    return P2

def ec_has_zero_coordinate(P: ECPoint) -> bool:
    return fp2_is_zero(P.x) or fp2_is_zero(P.z)




def ec_ladder3pt(
    m: int,
    P: ECPoint,
    Q: ECPoint,
    PQ: ECPoint,
    E: ECCurve
) -> tuple[bool, ECPoint]:
    """
    Functional equivalent of C ec_ladder3pt.

    Returns:
      (ok, R)
      ok = True on success, False on failure
      R  = P + m*Q (projective x-only), valid only if ok is True

    Notes:
      - C iterates NWORDS_ORDER * RADIX bits (e.g., 4*64 = 256 bits).
      - Here m is a Python int. We iterate either:
          * kbits bits if provided, else
          * NWORDS_ORDER*RADIX bits by default.
    """

    assert E.is_A24_computed_and_normalized

    if not fp2_is_one(E.A24.z):
        return False, P  # dummy

    # formulas invalid if PQ has zero coordinate
    if ec_has_zero_coordinate(PQ):
        return False, P  # dummy

    X0 = ec_point_init()
    X0 = copy_ec_point(Q)
    X1 = ec_point_init()
    X1 = copy_ec_point(P)
    X2 = ec_point_init()
    X2 = copy_ec_point(PQ)


    # EXACT C LOOP: word-by-word, LSB-first
    for i in range(NWORDS_ORDER):
        word = (m >> (i * RADIX)) & ((1 << RADIX) - 1)
        t = 1

        for _ in range(RADIX):
            bit_is_one = (word & t) != 0

            # C: swap if bit == 0
            swap = 0 if bit_is_one else 1

            X1, X2 = cswap_points(X1, X2, swap)
            X0, X1 = xDBLADD(X0, X1, X2, E.A24, True)
            X1, X2 = cswap_points(X1, X2, swap)

            t <<= 1

    return True, X1




def ec_normalize_point(
    P: ECPoint
) -> ECPoint:
    """
    Functional equivalent of:
        void ec_normalize_point(ec_point_t *P)
    """

    z_inv = fp2_inv(P.z)
    x = fp2_mul(P.x, z_inv)

    return ECPoint(
        x=x,
        z=fp2_set_one(),
    )

def AC_to_A24(curve: ECCurve) -> ECPoint:
    """
    Functional equivalent of:
        void AC_to_A24(ec_point_t *A24, const ec_curve_t *E)

    Returns a new ECPoint A24 = (A + 2C : 4C)
    """

    # If already computed, just return it
    if curve.is_A24_computed_and_normalized:
        return curve.A24

    # z = 2C
    z = fp2_add(curve.C, curve.C)

    # x = A + 2C
    x = fp2_add(curve.A, z)

    # z = 4C
    z = fp2_add(z, z)

    return ECPoint(x=x, z=z)


def ec_curve_normalize_A24(
    curve: ECCurve
) -> ECCurve:
    """
    Functional equivalent of:
        void ec_curve_normalize_A24(ec_curve_t *E)
    """

    if curve.is_A24_computed_and_normalized:
        return curve
    A24 = AC_to_A24(curve)
    A24 = ec_normalize_point(A24)

    assert fp2_is_one(A24.z)

    return ECCurve(
        A=curve.A,
        C=curve.C,
        A24=A24,
        is_A24_computed_and_normalized=True,
    )

def cswap_points(P: ECPoint, Q: ECPoint, option: int) -> tuple[ECPoint, ECPoint]:
    """
    Functional equivalent of constant-time cswap_points.

    option == 0  -> return (P, Q)
    option != 0  -> return (Q, P)
    """

    if option:
        return Q, P
    else:
        return P, Q

def xDBLADD(
    P: ECPoint,
    Q: ECPoint,
    PQ: ECPoint,
    A24: ECPoint,
    A24_normalized: bool
) -> tuple[ECPoint, ECPoint]:
    """
    Functional equivalent of C xDBLADD.

    Inputs:
      P, Q: ECPoint  (projective x-only)
      PQ:   ECPoint  (difference P-Q)
      A24:  ECPoint  (A24 = (A+2C : 4C) or (A+2C/4C : 1) if normalized)
      A24_normalized: bool

    Outputs:
      (R, S) where:
        R = 2*P
        S = P+Q
    """

    # t0 = XP + ZP
    t0 = fp2_add(P.x, P.z)

    # t1 = XP - ZP
    t1 = fp2_sub(P.x, P.z)

    # Rx = t0^2
    Rx = fp2_sqr(t0)

    # t2 = XQ - ZQ
    t2 = fp2_sub(Q.x, Q.z)

    # Sx = XQ + ZQ
    Sx = fp2_add(Q.x, Q.z)

    # t0 = t0 * t2
    t0 = fp2_mul(t0, t2)

    # Rz = t1^2
    Rz = fp2_sqr(t1)

    # t1 = t1 * Sx
    t1 = fp2_mul(t1, Sx)

    # t2 = Rx - Rz
    t2 = fp2_sub(Rx, Rz) # op9

    # if not normalized: Rz = Rz * A24.z
    # if not A24_normalized:
    #     Rz = fp2_mul(Rz, A24.z)

    # Rx = Rx * Rz
    Rx = fp2_mul(Rx, Rz) # op10

    # Sx = A24.x * t2
    Sx2 = fp2_mul(A24.x, t2) # op11

    # Sz = t0 - t1
    Sz = fp2_sub(t0, t1) # op12

    # Rz = Rz + Sx
    Rz = fp2_add(Rz, Sx2) # op13

    # Sx = t0 + t1
    Sx = fp2_add(t0, t1) # op14

    # Rz = Rz * t2
    Rz = fp2_mul(Rz, t2) # op15

    # Sz = Sz^2
    Sz = fp2_sqr(Sz) # op16

    # Sx = Sx^2
    Sx = fp2_sqr(Sx) # op17

    # Sz = Sz * PQ.x
    Sz = fp2_mul(Sz, PQ.x)

    # Sx = Sx * PQ.z
    Sx = fp2_mul(Sx, PQ.z)

    R = ECPoint(x=Rx, z=Rz)
    S = ECPoint(x=Sx, z=Sz)
    return R, S


def xMUL(P: ECPoint, k: int, kbits: int, curve: ECCurve) -> ECPoint:
    """
    Functional equivalent of C xMUL (Montgomery ladder).

    Inputs:
      P: ECPoint (projective x-only)
      scalar: int (Python integer scalar)
      kbits: int (bitlength to iterate)
      curve: ECCurve

    Output:
      ECPoint Q = scalar * P (projective x-only)
    """

    # Build A24
    # if not curve.is_A24_computed_and_normalized:
    #     # A24 = (A + 2C : 4C)
    #     print("not normalized")
    #     z = fp2_add(curve.C, curve.C)      # 2C
    #     z4 = fp2_add(z, z)                 # 4C
    #     x = fp2_add(z, curve.A)            # A + 2C
    #     A24 = ECPoint(x=x, z=z4)
    # else:
    A24 = curve.A24
        # assert fp2_is_one(A24.z)

    # R0 = (1:0), R1 = P
    R0 = ec_point_init()       # returns identity point (1:0)
    R1 = ec_point_init() 
    R1 = copy_ec_point(P)         # functional copy

    prevbit = 0

    # Main ladder loop
    for i in range(kbits - 1, -1, -1):
        bit = (k >> i) & 1
        swap = bit ^ prevbit
        prevbit = bit

        # Oracle swap (C uses mask = 0 - swap)
        R0, R1 = cswap_points(R0, R1, swap)

        # Differential double-and-add
        R0, R1 = xDBLADD(R0, R1, P, A24, True)


    # Final swap (swap = 0 ^ prevbit)
    swap = prevbit
    R0, R1 = cswap_points(R0, R1, swap)

    return R0




def ec_mul(
    scalar: int,
    kbits: int,
    P: ECPoint,
    curve: ECCurve
) -> ECPoint:
    """
    Functional equivalent of:
        void ec_mul(ec_point_t *res, const digit_t *scalar,
                    int kbits, const ec_point_t *P, ec_curve_t *curve)
    """
    if kbits > 50:
        curve = ec_curve_normalize_A24(curve)

    return xMUL(P, scalar, kbits, curve)

def xDBL(
    P: ECPoint,
    curve: ECCurve
) -> ECPoint:
    """
    Functional equivalent of C xDBL.
    Uses (A:C) directly.
    """

    t0 = fp2_add(P.x, P.z)
    t0 = fp2_sqr(t0)

    t1 = fp2_sub(P.x, P.z)
    t1 = fp2_sqr(t1)

    t2 = fp2_sub(t0, t1)

    t3 = fp2_add(curve.C, curve.C)      # 2C
    t1 = fp2_mul(t1, t3)
    t1 = fp2_add(t1, t1)                # 4C * (XP-ZP)^2

    Rx = fp2_mul(t0, t1)

    t0 = fp2_add(t3, curve.A)           # A + 2C
    t0 = fp2_mul(t0, t2)
    t0 = fp2_add(t0, t1)

    Rz = fp2_mul(t0, t2)

    R = ECPoint(x=Rx, z=Rz)

    return R




def xDBL_A24(
    P: ECPoint,
    A24: ECPoint,
    A24_normalized: bool
) -> ECPoint:
    """
    Functional equivalent of C xDBL_A24.
    """

    t0 = fp2_add(P.x, P.z)
    t0 = fp2_sqr(t0)

    t1 = fp2_sub(P.x, P.z)
    t1 = fp2_sqr(t1)

    t2 = fp2_sub(t0, t1)

    if not A24_normalized:
        t1 = fp2_mul(t1, A24.z)


    X = fp2_mul(t0, t1)

    t0 = fp2_mul(t2, A24.x)
    

    t0 = fp2_add(t0, t1)

    Z = fp2_mul(t0, t2)

    return ECPoint(x=X, z=Z)


def ec_dbl( P: ECPoint, curve: ECCurve):

    res = ec_point_init()
    # If A24 = ((A+2)/4 : 1) we save multiplications
    if curve.is_A24_computed_and_normalized:
        assert fp2_is_one(curve.A24.z)
        res = xDBL_A24(P, curve.A24, True)
    else:
        #Otherwise we compute A24 on the fly for doubling
        res = xDBL(P, curve)
    
    return res

def ec_dbl_iter(res: ECPoint, n: int, P: ECPoint, curve: ECCurve):
    """
    Perform n successive elliptic curve doublings starting from P.

    Parameters:
        res   : ec_point_t (output, modified in-place)
        n     : int (number of doublings)
        P     : ec_point_t (input point)
        curve : ec_curve_t
    """

    curve_new = curve
    if n == 0:
        res = copy_ec_point(P)
        return res, curve_new

    # When the chain is long enough, normalize A24
    if n > 50:
        curve_new = ec_curve_normalize_A24(curve_new)
    
    # If A24 is already computed and normalized, use the optimized doubling
    if curve_new.is_A24_computed_and_normalized:
        assert fp2_is_one(curve_new.A24.z)

        # First doubling
        res = xDBL_A24(P, curve_new.A24, True)

        # Remaining doublings
        for _ in range(n - 1):
            assert fp2_is_one(curve_new.A24.z)
            res = xDBL_A24(res, curve_new.A24, True)

    else:
        # Fallback: generic doubling
        res = xDBL(P, curve_new)

        for _ in range(n - 1):
            res = xDBL(res, curve_new)

    return res, curve_new

def ec_dbl_iter_basis(res: ECBasis, n: int, B: ECBasis, curve: ECCurve):
    res_new = res
    new_curve = curve
    res_new.P, new_curve = ec_dbl_iter(res_new.P, n, B.P, curve)
    res_new.Q, new_curve = ec_dbl_iter(res_new.Q, n, B.Q, curve)
    res_new.PmQ, new_curve = ec_dbl_iter(res_new.PmQ, n, B.PmQ, curve)

    return res_new, new_curve


def ec_is_zero(P: ECPoint) -> bool:
    """
    Functional equivalent of C ec_is_zero.
    Returns True iff P is the point at infinity.
    """
    return fp2_is_zero(P.z)

def ec_is_equal(P: ECPoint, Q: ECPoint):
    # Evaluate if two points in Montgomery coordinates (X:Z) are equal
    # Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    t0 = Fp2()
    t1 = Fp2()

    # Check if P, Q are the points at infinity
    l_zero = ec_is_zero(P)
    r_zero = ec_is_zero(Q)

    # Check if PX * QZ = QX * PZ
    t0 = fp2_mul(P.x, Q.z)
    t1 = fp2_mul(P.z, Q.x)
    lr_equal = fp2_is_equal(t0, t1)

    # Points are equal if
    # - Both are zero, or
    # - neither are zero AND PX * QZ = QX * PZ
    return (l_zero and r_zero) or ((not l_zero) & (not r_zero) * lr_equal)


def ec_is_two_torsion(P: ECPoint, E: ECCurve) -> bool:
    """
    Functional equivalent of C ec_is_two_torsion.
    Returns True if P is a 2-torsion point.
    """

    if ec_is_zero(P):
        return False

    t0 = fp2_add(P.x, P.z)
    t0 = fp2_sqr(t0)

    t1 = fp2_sub(P.x, P.z)
    t1 = fp2_sqr(t1)

    t2 = fp2_sub(t0, t1)

    t1 = fp2_add(t0, t1)

    t2 = fp2_mul(t2, E.A)
    t1 = fp2_mul(t1, E.C)
    t1 = fp2_add(t1, t1)

    # t0 = 4 * (C*X^2 + C*Z^2 + A*X*Z)
    t0 = fp2_add(t1, t2)

    x_is_zero = fp2_is_zero(P.x)
    tmp_is_zero = fp2_is_zero(t0)

    # two torsion if x == 0 or x^2 + A*x + 1 == 0
    return x_is_zero or tmp_is_zero

def ec_is_four_torsion(P: ECPoint, E: ECCurve) -> bool:
    """
    Functional equivalent of C ec_is_four_torsion.
    """

    test = xDBL_A24(P, E.A24, E.is_A24_computed_and_normalized)
    return ec_is_two_torsion(test, E)

def mp_shiftr(a: int, shift: int):
    return a & 1, a>>1


def xADD(P: ECPoint, Q: ECPoint, PQ: ECPoint) -> ECPoint:
    """
    Differential addition of Montgomery points in projective coordinates (X:Z).

    Input:
      P  = (XP : ZP)
      Q  = (XQ : ZQ)
      PQ = P - Q = (XPQ : ZPQ)

    Output:
      R = P + Q = (XR : ZR)
    """

    t0 = fp2_add(P.x, P.z)
    t1 = fp2_sub(P.x, P.z)
    t2 = fp2_add(Q.x, Q.z)
    t3 = fp2_sub(Q.x, Q.z)

    t0 = fp2_mul(t0, t3)
    t1 = fp2_mul(t1, t2)

    t2 = fp2_add(t0, t1)
    t3 = fp2_sub(t0, t1)

    t2 = fp2_sqr(t2)
    t3 = fp2_sqr(t3)

    Rx = fp2_mul(PQ.z, t2)
    Rz = fp2_mul(PQ.x, t3)

    R = ECPoint(x=Rx, z=Rz)

    return R



def xDBLMUL( P: ECPoint, k: int, Q: ECPoint, l: int, PQ: ECPoint, kbits: int, curve: ECCurve):
    S = ec_point_init()

    if ec_has_zero_coordinate(P) or ec_has_zero_coordinate(Q) or ec_has_zero_coordinate(PQ):
        return 0

    r = [0] * (2 * kbits)

    bitk0 = (k & 1)
    bitl0 = (l & 1)

    maskk = 0 - bitk0
    maskl = 0 - bitl0

    sigma0 = (bitk0 ^ 1)
    sigma1 = (bitl0 ^ 1)

    evens = sigma0 + sigma1
    mevens = 0 - (evens & 1)

    sigma0 = (sigma0 & mevens)
    sigma1 = (sigma1 & mevens) | (1 & (~mevens))

    one0 = 1

    k_t = k - one0
    l_t = l - one0

    if maskk == 0:
        k_t = k_t
    else:
        k_t = k

    if maskl == 0:
        l_t = l_t
    else:
        l_t = l


    pre_sigma = 0

    for i in range(kbits):
        maskk = 0 - (sigma0 ^ pre_sigma)
        if maskk == 0:
            k_t = k_t
            l_t = l_t
        else:
            temp = k_t
            k_t = l_t
            l_t = temp

        if i == kbits - 1:
            bs1_ip1 = 0
            bs2_ip1 = 0
        else:
            bs1_ip1, k_t = mp_shiftr(k_t, 1)
            bs2_ip1, l_t = mp_shiftr(l_t, 1)
        

        bs1_i = k_t & 1
        bs2_i = l_t & 1

        r[2*i]   = bs1_i ^ bs1_ip1
        r[2*i+1] = bs2_i ^ bs2_ip1

        pre_sigma = sigma0
        maskk = 0 - r[2 * i + 1]
        if maskk == 0:
            sigma0, sigma1 = sigma0, sigma1
        else:
            sigma0, sigma1 = sigma1, sigma0

    R0 = ec_point_init()
    maskk = 0 - sigma0
    if maskk == 0:        
        R1 = P
        R2 = Q
    else:
        R1 = Q
        R2 = P

    DIFF1a = ec_point_init()
    DIFF1b = ec_point_init()

    DIFF1a = copy_ec_point(R1)
    DIFF1b = copy_ec_point(R2)

    R2 = xADD(R1, R2, PQ)
    if ec_has_zero_coordinate(R2):
        return 0

    DIFF2a = copy_ec_point(R2)
    DIFF2b = copy_ec_point(PQ)

    A_is_zero = fp2_is_zero(curve.A)

    for i in reversed(range(kbits)):
        h = r[2*i] + r[2*i+1]

        maskk = 0 - (h & 1)
        if maskk == 0:
            T0 = R0
        else:
            T0 = R1

        maskk = 0 - (h >> 1)
        if maskk == 0:
            T0 = T0
        else:
            T0 = R2

        if A_is_zero:
            T0 = xDBL_E0(T0)
        else:
            assert fp2_is_one(curve.A24.z)
            T0 = xDBL_A24(T0, curve.A24, True)
        
        maskk = 0 - r[2 * i + 1]; 
        if maskk == 0:
            T1 = R0
            T2 = R1
        else:
            T1 = R1
            T2 = R2

        DIFF1a, DIFF1b = cswap_points(DIFF1a, DIFF1b, maskk)

        T1 = xADD(T1, T2, DIFF1a)
        T2 = xADD(R0, R2, DIFF2a)

        maskk = 0 - (h & 1)
        DIFF2a, DIFF2b = cswap_points(DIFF2a, DIFF2b, maskk)

        R0 = copy_ec_point(T0)
        R1 = copy_ec_point(T1)
        R2 = copy_ec_point(T2)


    if mevens == 0:
        S = R0
    else:
        S = R1
    
    maskk = 0 - (bitk0 & bitl0)
    if maskk == 0:
        S = S
    else:
        S = R2

    return S, 1




def ec_biscalar_mul(scalarP: int,
                    scalarQ: int,
                    kbits: int,
                    PQ: ECBasis,
                    curve: ECCurve) -> int:

    if fp2_is_zero(PQ.PmQ.z):
        return 0

    # Special case: kbits == 1
    if kbits == 1:
        if (not ec_is_two_torsion(PQ.P, curve) or
            not ec_is_two_torsion(PQ.Q, curve) or
            not ec_is_two_torsion(PQ.PmQ, curve)):
            return 0

        bP = scalarP & 1
        bQ = scalarQ & 1

        if bP == 0 and bQ == 0:
            res = ec_point_init()        # (1:0)
        elif bP == 1 and bQ == 0:
            res = copy_ec_point(PQ.P)
        elif bP == 0 and bQ == 1:
            res = copy_ec_point(PQ.Q)
        elif bP == 1 and bQ == 1:
            res = copy_ec_point(PQ.PmQ)
        else:
            raise AssertionError("Impossible scalar combination")

        return 1

    # General case
    E = ECCurve(
        A=fp2_copy(curve.A),
        C=fp2_copy(curve.C),
        A24=curve.A24,
        is_A24_computed_and_normalized=curve.is_A24_computed_and_normalized
    )

    if not fp2_is_zero(curve.A):
        E = ec_curve_normalize_A24(E)

    res, ret = xDBLMUL(
        PQ.P, scalarP,
        PQ.Q, scalarQ,
        PQ.PmQ,
        kbits,
        E
    )
    return res, ret

def jac_to_ws(ao3: Fp2, P: jac_point_t, curve: ECCurve)->any:
    one = ONE

    ao3_new = ao3

    Q = jac_point_t()

    if not fp2_is_zero(curve.A):
        ao3_new.re = fp_div3( curve.A.re, q, True)
        ao3_new.im = fp_div3( curve.A.im, q, True)

        t = fp2_sqr(P.z)
        Q.x = fp2_mul(ao3_new, t)
        Q.x = fp2_add( Q.x, P.x)

        t = fp2_sqr(t)
        a = fp2_set_zero()
        a = fp2_mul( ao3_new, curve.A)

        a.re = fp_sub(one, a.re)
        a.im = fp_neg(a.im)

        t = fp2_mul( t, a)
    else:
        Q.x = fp2_copy(P.x)
        t = fp2_sqr(P.z)
        t = fp2_sqr( t)

    Q.y = fp2_copy(P.y)
    Q.z = fp2_copy(P.z)

    return Q, t, ao3_new


def jac_from_ws(Q: jac_point_t, P: jac_point_t, ao3: Fp2, curve: ECCurve):
    Q_new = jac_point_t()
    Q_new.x = Q.x
    Q_new.y = Q.y
    Q_new.z = Q.z
    if not fp2_is_zero(curve.A):
        t   = fp2_set_zero()
        t   = fp2_sqr(P.z)
        t   = fp2_mul(t, ao3)
        Q_new.x = fp2_sub(P.x, t)

    Q_new.y = fp2_copy(P.y)
    Q_new.z = fp2_copy(P.z)

    return Q_new



def DBL(P: jac_point_t, AC: ECPoint):
    Q = jac_point_t()
    t0 = fp2_set_zero()
    t1 = fp2_set_zero()
    t2 = fp2_set_zero()
    t3 = fp2_set_zero()

    flag = fp2_is_zero(P.x) & fp2_is_zero(P.z)

    t0 = fp2_sqr(P.x)
    t1 = fp2_add(t0, t0)
    t0 = fp2_add(t0, t1)

    t1 = fp2_sqr( P.z)
    t2 = fp2_mul( P.x, AC.A)
    t2 = fp2_add( t2, t2)
    t2 = fp2_add(t1, t2)
    t2 = fp2_mul(t1, t2)
    t2 = fp2_add(t0, t2)

    Q.z = fp2_mul(P.y, P.z)
    Q.z = fp2_add(Q.z, Q.z)

    t0 = fp2_sqr(Q.z)
    t0 = fp2_mul(t0, AC.A)

    t1 = fp2_sqr(P.y)
    t1 = fp2_add(t1, t1)

    t3 = fp2_add(P.x, P.x)
    t3 = fp2_mul(t1, t3)

    Q.x = fp2_sqr(t2)
    Q.x = fp2_sub(Q.x, t0)
    Q.x = fp2_sub(Q.x, t3)
    Q.x = fp2_sub(Q.x, t3)

    Q.y = fp2_sub(t3, Q.x)
    Q.y = fp2_mul(Q.y, t2)

    t1 = fp2_sqr(t1)
    Q.y = fp2_sub( Q.y, t1)
    Q.y = fp2_sub( Q.y, t1)

    ctl = -flag & 0xFFFFFFFF
    Q.z = fp2_select(Q.z, P.z, ctl)
    Q.x = fp2_select(Q.x, P.x, ctl)

    return Q





def double_couple_jac_point(inp: theta_couple_jac_point_t, E1E2: theta_couple_curve_t):
    out = theta_couple_jac_point_t(P1=None, P2=None)
    out.P1 = DBL(inp.P1, E1E2.E1)
    out.P2 = DBL(inp.P2, E1E2.E2)

    return out

def DBLW(P: jac_point_t, t: Fp2):
    Q = jac_point_t()
    u = Fp2()
    flag = fp2_is_zero(P.x) & fp2_is_zero(P.z)

    xx = fp2_set_zero()
    c = fp2_set_zero()
    cc = fp2_set_zero()
    r = fp2_set_zero()
    s = fp2_set_zero()
    m = fp2_set_zero()

    xx = fp2_sqr(P.x)

    c = fp2_sqr( P.y)
    c = fp2_add( c, c)

    cc  = fp2_sqr( c)
    r   = fp2_add(cc, cc)

    s = fp2_add(P.x, c)
    s = fp2_sqr(s)
    s = fp2_sub(s, xx)
    s = fp2_sub(s, cc)

    m = fp2_add(xx, xx)
    m = fp2_add(m, xx)
    m = fp2_add(m, t)

    Q.x = fp2_sqr(m)
    Q.x = fp2_sub(Q.x, s)
    Q.x = fp2_sub(Q.x, s)

    Q.z = fp2_mul(P.y, P.z)
    Q.z = fp2_add(Q.z, Q.z)

    Q.y = fp2_sub(s, Q.x)
    Q.y = fp2_mul(Q.y, m)
    Q.y = fp2_sub(Q.y, r)

    u = fp2_mul(t, r)
    u = fp2_add(u, u)

    ctl = -flag & 0xFFFFFFFF
    Q.x = fp2_select(Q.x, P.x, ctl)
    Q.z = fp2_select(Q.z, P.z, ctl)

    return Q, u



def double_couple_jac_point_iter(n: int, inp: theta_couple_jac_point_t, E1E2: theta_couple_curve_t):
    out = theta_couple_jac_point_t(P1=jac_point_t(), P2=jac_point_t())
    if n == 0:
        out.P1 = copy_jac_point(inp.P1)
        out.P2 = copy_jac_point(inp.P2)
    elif n == 1:
        out = double_couple_jac_point(inp, E1E2)
    else:
        t1 = fp2_set_zero()
        t2 = fp2_set_zero()
        a1 = fp2_set_zero()
        a2 = fp2_set_zero()

        out.P1, t1, a1 = jac_to_ws(a1, inp.P1, E1E2.E1)
        out.P2, t2, a2 = jac_to_ws(a2, inp.P2, E1E2.E2)

        out.P1, t1 = DBLW( out.P1, t1)
        out.P2, t2 = DBLW( out.P2, t2)

        for _ in range(n - 1):
            out.P1, t1 = DBLW(out.P1, t1)
            out.P2, t2 = DBLW(out.P2, t2)

        out.P1 = jac_from_ws(out.P1, out.P1, a1, E1E2.E1)
        out.P2 = jac_from_ws(out.P2, out.P2, a2, E1E2.E2)


    return out



