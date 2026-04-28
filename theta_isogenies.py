from arith_ops import *
from ec import *
from basis import *

ENABLE_SIGN = True


def print_hex_point(name, point):
    print("---" + str(name) + "---")
    print("A.re: ", hex(point.x.re))
    print("A.im: ", hex(point.x.im))
    print("Z.re: ", hex(point.z.re))
    print("Z.im: ", hex(point.z.im))


def print_fp2(name: str, A : Fp2):
    print("---" + str(name) + "---")
    print("A.re: ", hex(A.re))
    print("A.im: ", hex(A.im))

def print_curve(name: str ,curve: ECCurve):
    print("---" + str(name) + "---")
    print_fp2("curve.A: ",  curve.A)
    print_fp2("curve.C: ",  curve.C)
    if curve.A24 != None:
        print_hex_point("curve.A24: ", curve.A24)
    else:
        print("curve.A24: ", None)
    print("curve.normalized: ", curve.is_A24_computed_and_normalized)


def fp2_sqrt_verify(a: Fp2):

    tmp = a
    t0 = fp2_copy(a)
    t0 = fp2_sqrt(t0)
    t1 = fp2_sqr(t0)

    return t0, fp2_is_equal(a, t1)


def ec_recover_y(Px: Fp2, curve: ECCurve):
    #Recover y-coordinate of a point on the Montgomery curve y^2 = x^3 + Ax^2 + x
    t0 = fp2_sqr(Px)
    y = fp2_mul(t0, curve.A); #  Ax^2
    y = fp2_add(y, Px);      #  Ax^2 + x
    t0 = fp2_mul(t0, Px)
    y = fp2_add(y, t0) #  x^3 + Ax^2 + x
    # This is required, because we do not yet know that our curves are
    # supersingular so our points live on the twist with B = 1.
    y, ret = fp2_sqrt_verify(y)
    return y , ret


def lift_basis_normalized(B: ECBasis,
                          E: ECCurve) -> int:
    # Preconditions
    assert fp2_is_one(B.P.z)
    assert fp2_is_one(E.C)

    P = jac_point_t(x=Fp2(re=0, im=0), y=Fp2(re=0, im=0),z=Fp2(re=0, im=0))
    Q = jac_point_t(x=Fp2(re=0, im=0), y=Fp2(re=0, im=0),z=Fp2(re=0, im=0))

    # Copy x-coordinates
    P.x = fp2_copy(B.P.x)
    Q.x = fp2_copy(B.Q.x)
    Q.z = fp2_copy(B.Q.z)

    P.z = fp2_set_one()

    # Recover y-coordinate of P
    P.y, ret = ec_recover_y(P.x, E)

    # --- Okeya–Sakurai algorithm to recover y(Q) ---

    v1 = fp2_mul(P.x, Q.z)
    v2 = fp2_add(Q.x, v1)
    v3 = fp2_sub(Q.x, v1)

    v3 = fp2_sqr(v3)
    v3 = fp2_mul(v3, B.PmQ.x)

    v1 = fp2_add(E.A, E.A)       # 2A
    v1 = fp2_mul(v1, Q.z)
    v2 = fp2_add(v2, v1)

    v4 = fp2_mul(P.x, Q.x)
    v4 = fp2_add(v4, Q.z)
    v2 = fp2_mul(v2, v4)

    v1 = fp2_mul(v1, Q.z)
    v2 = fp2_sub(v2, v1)
    v2 = fp2_mul(v2, B.PmQ.z)

    Q.y = fp2_sub(v3, v2)

    v1 = fp2_add(P.y, P.y)
    v1 = fp2_mul(v1, Q.z)
    v1 = fp2_mul(v1, B.PmQ.z)

    Q.x = fp2_mul(Q.x, v1)
    Q.z = fp2_mul(Q.z, v1)

    # Convert to Jacobian coordinates
    v1 = fp2_sqr(Q.z)
    Q.y = fp2_mul(Q.y, v1)
    Q.x = fp2_mul(Q.x, Q.z)

    return P, Q, ret


def lift_basis(B: ECBasis,
               E: ECCurve) -> int:
    # Normalize curve and point:
    #   (A : C) -> (A/C : 1)
    #   (X : Z) -> (X/Z : 1)

    P = jac_point_t()
    Q = jac_point_t()

    inverses = [
        fp2_copy(B.P.z),
        fp2_copy(E.C),
    ]


    # assert inverses[0].re in dict_lift, "ERROR OF FP2_INV" 
    # assert inverses[1].re in dict_lift, "ERROR OF FP2_INV" 
    # inverses[0] = dict_lift[inverses[0].re]
    # inverses[1] = dict_lift[inverses[1].re]
    inverses[0] = fp2_inv(inverses[0])
    inverses[1] = fp2_inv(inverses[1])
    #inverses = fp2_batched_inv(inverses)

    B.P.z = fp2_set_one()
    E.C   = fp2_set_one()

    B.P.x = fp2_mul(B.P.x, inverses[0])
    E.A   = fp2_mul(E.A, inverses[1])

    # Lift the basis
    P, Q, ret = lift_basis_normalized(B, E)
    
    return P, Q, ret


def jac_to_xz(xyP: jac_point_t):
    P = ec_point_init()
    P.x = fp2_copy(xyP.x)
    P.z = fp2_copy(xyP.z)
    P.z = fp2_sqr(P.z)

    # If xyP = (0:1:0), we currently have P=(0 : 0) but we want to set P=(1:0)
    one = fp2_set_one()

    c1 = fp2_is_zero(P.x)
    c2 = fp2_is_zero(P.z)
    P.x = fp2_select(P.x, one, c1 and c2)

    return P



def couple_jac_to_xz(xyP: theta_couple_jac_point_t):
    P = theta_couple_point_t()
    P.P1 = jac_to_xz(xyP.P1)
    P.P2 = jac_to_xz(xyP.P2)

    return P


def action_by_translation_z_and_det(P4: ECPoint, 
                                    P2: ECPoint):
    """
    Collect Z coordinate and determinant for batched inversion.
    """

    z_inv   = Fp2()
    det_inv = Fp2()

    # z_inv = P4.z
    z_inv = fp2_copy(P4.z)   # or: z_inv = P4.z.copy()

    # det = x4 * z2 - z4 * x2
    det = fp2_mul(P4.x, P2.z)
    tmp = fp2_mul(P4.z, P2.x)
    det = fp2_sub(det, tmp)

    det_inv = fp2_copy(det)

    return z_inv, det_inv



def double_couple_point(in1: theta_couple_point_t, E1E2: theta_couple_curve_t):
    out =  theta_couple_point_t()
    out.P1 = ec_dbl(in1.P1, E1E2.E1)
    out.P2 = ec_dbl(in1.P2, E1E2.E2)

    return out

def verify_two_torsion(K1_2: theta_couple_point_t, 
                       K2_2 :theta_couple_point_t, 
                       E12  :theta_couple_curve_t):

    # First check if any point in K1_2 or K2_2 is zero, if they are then the points did not have
    # order 8 when we started gluing
    if (ec_is_zero(K1_2.P1) or ec_is_zero(K1_2.P2) or ec_is_zero(K2_2.P1) or ec_is_zero(K2_2.P2)):
        return False

    # Now ensure that P1, Q1 and P2, Q2 are independent. For points of order two this means
    # that they're not the same
    if ec_is_equal(K1_2.P1, K2_2.P1) or ec_is_equal(K1_2.P2, K2_2.P2):
        return False
    

    # Finally, double points to ensure all points have order exactly 0
    O1 = theta_couple_point_t()
    O2 = theta_couple_point_t()
    O1 = double_couple_point(K1_2, E12)
    O2 = double_couple_point(K2_2, E12)
    # If this check fails then the points had order 2*f for some f, and the kernel is malformed.
    if not (ec_is_zero(O1.P1) and ec_is_zero(O1.P2) and ec_is_zero(O2.P1) and ec_is_zero(O2.P2)):
        return False

    return True



def action_by_translation_compute_matrix(P4: ECPoint,
                                         P2: ECPoint,
                                         z_inv: Fp2,
                                         det_inv: Fp2):

    tmp = Fp2()

    G = translation_matrix_t()

    # Gi.g10 = uij xij /detij - xij/zij
    tmp     = fp2_mul(P4.x, z_inv)
    G.g10   = fp2_mul(P4.x, P2.x)
    G.g10   = fp2_mul(G.g10, det_inv)
    G.g10   = fp2_sub(G.g10, tmp)

    # Gi.g11 = uij zij * detij
    G.g11 = fp2_mul(P2.x, det_inv)
    G.g11 = fp2_mul(G.g11, P4.z)

    #3Gi.g00 = -Gi.g11
    G.g00 = fp2_neg(G.g11)

    # Gi.g01 = - wij zij detij
    G.g01 = fp2_mul(P2.z, det_inv)
    G.g01 = fp2_mul(G.g01, P4.z)
    G.g01 = fp2_neg(G.g01)

    return G



def action_by_translation(K1_4: theta_couple_point_t, 
                          K2_4: theta_couple_point_t, 
                          E12: theta_couple_curve_t):
    """
    Gi: list[translation_matrix_t] of length 4
    """
    Gi = [translation_matrix_t() for i in range(4)]

    # Compute points of order 2 from Ki_4
    K1_2 = theta_couple_point_t(P1=jac_point_t(), P2=jac_point_t())
    K2_2 = theta_couple_point_t(P1=jac_point_t(), P2=jac_point_t())

    K1_2 = double_couple_point(K1_4, E12)
    K2_2 = double_couple_point(K2_4, E12)

    ret = verify_two_torsion(K1_2, K2_2, E12)

    if not ret:
        return False

    # We invert four Z coordinates and four determinants (batched)
    inverses = [fp2_set_zero() for _ in range(8)]

    inverses[0], inverses[4] = action_by_translation_z_and_det(K1_4.P1, K1_2.P1)
    inverses[1], inverses[5] = action_by_translation_z_and_det(K1_4.P2, K1_2.P2)
    inverses[2], inverses[6] = action_by_translation_z_and_det(K2_4.P1, K2_2.P1)
    inverses[3], inverses[7] = action_by_translation_z_and_det(K2_4.P2, K2_2.P2)

    #inverses = fp2_batched_inv(inverses, 8)
    inverses = [fp2_inv(inverses[i]) for i in range(len(inverses))]

    if fp2_is_zero(inverses[0]):
        return False  # invalid input (matches C behavior)

    Gi[0] = action_by_translation_compute_matrix(K1_4.P1, K1_2.P1, inverses[0], inverses[4])
    Gi[1] = action_by_translation_compute_matrix(K1_4.P2, K1_2.P2, inverses[1], inverses[5])
    Gi[2] = action_by_translation_compute_matrix(K2_4.P1, K2_2.P1, inverses[2], inverses[6])
    Gi[3] = action_by_translation_compute_matrix(K2_4.P2, K2_2.P2, inverses[3], inverses[7])

    return True, Gi



def gluing_change_of_basis( K1_4: theta_couple_point_t,
                            K2_4: theta_couple_point_t,
                            E12 : theta_couple_curve_t ):

    M = basis_change_matrix_init()
    # Compute the four 2x2 matrices for the action by translation
    # on the four points:
    Gi = [translation_matrix_t() for i in range(4)]
    ret, Gi = action_by_translation(K1_4, K2_4, E12)
    if not ret:
        return 0

    # Computation of the 4x4 matrix from Mij
    # t001, t101 (resp t002, t102) first column of M11 * M21 (resp M12 * M22)
    # fp2_t t001, t101, t002, t102, tmp

    t001 = fp2_mul(Gi[0].g00, Gi[2].g00)
    tmp  = fp2_mul(Gi[0].g01, Gi[2].g10)
    t001 = fp2_add(t001, tmp)

    t101 = fp2_mul(Gi[0].g10, Gi[2].g00)
    tmp = fp2_mul(Gi[0].g11, Gi[2].g10)
    t101 = fp2_add(t101, tmp)

    t002 = fp2_mul(Gi[1].g00, Gi[3].g00)
    tmp  = fp2_mul(Gi[1].g01, Gi[3].g10)
    t002 = fp2_add(t002, tmp)

    t102 = fp2_mul(Gi[1].g10, Gi[3].g00)
    tmp  = fp2_mul(Gi[1].g11, Gi[3].g10)
    t102 = fp2_add(t102, tmp)

    # trace for the first row
    M.m[0][0] = fp2_set_one()
    M.m[0][0] = fp2_add(M.m[0][0], fp2_mul(t001, t002))
    M.m[0][0] = fp2_add(M.m[0][0], fp2_mul(Gi[2].g00, Gi[3].g00))
    M.m[0][0] = fp2_add(M.m[0][0], fp2_mul(Gi[0].g00, Gi[1].g00))

    M.m[0][1] = fp2_mul(t001, t102)
    M.m[0][1] = fp2_add(M.m[0][1], fp2_mul(Gi[2].g00, Gi[3].g10))
    M.m[0][1] = fp2_add(M.m[0][1], fp2_mul(Gi[0].g00, Gi[1].g10))

    M.m[0][2] = fp2_mul(t101, t002)
    M.m[0][2] = fp2_add(M.m[0][2], fp2_mul(Gi[2].g10, Gi[3].g00))
    M.m[0][2] = fp2_add(M.m[0][2], fp2_mul(Gi[0].g10, Gi[1].g00))

    M.m[0][3] = fp2_mul(t101, t102)
    M.m[0][3] = fp2_add(M.m[0][3], fp2_mul(Gi[2].g10, Gi[3].g10))
    M.m[0][3] = fp2_add(M.m[0][3], fp2_mul(Gi[0].g10, Gi[1].g10))

    # Compute the action of (0,out.K2_4.P2) for the second row
    M.m[1][0] = fp2_mul(Gi[3].g00, M.m[0][0])
    M.m[1][0] = fp2_add(M.m[1][0], fp2_mul(Gi[3].g01, M.m[0][1]))

    M.m[1][1] = fp2_mul(Gi[3].g10, M.m[0][0])
    M.m[1][1] = fp2_add(M.m[1][1], fp2_mul(Gi[3].g11, M.m[0][1]))

    M.m[1][2] = fp2_mul(Gi[3].g00, M.m[0][2])
    M.m[1][2] = fp2_add(M.m[1][2], fp2_mul(Gi[3].g01, M.m[0][3]))

    M.m[1][3] = fp2_mul(Gi[3].g10, M.m[0][2])
    M.m[1][3] = fp2_add(M.m[1][3], fp2_mul(Gi[3].g11, M.m[0][3]))


    # compute the action of (K1_4.P1,0) for the third row
    M.m[2][0] = fp2_mul(Gi[0].g00, M.m[0][0])
    M.m[2][0] = fp2_add(M.m[2][0], fp2_mul(Gi[0].g01, M.m[0][2]))

    M.m[2][1] = fp2_mul(Gi[0].g00, M.m[0][1])
    M.m[2][1] = fp2_add(M.m[2][1], fp2_mul(Gi[0].g01, M.m[0][3]))

    M.m[2][2] = fp2_mul(Gi[0].g10, M.m[0][0])
    M.m[2][2] = fp2_add(M.m[2][2], fp2_mul(Gi[0].g11, M.m[0][2]))

    M.m[2][3] = fp2_mul(Gi[0].g10, M.m[0][1])
    M.m[2][3] = fp2_add(M.m[2][3], fp2_mul(Gi[0].g11, M.m[0][3]))


    # compute the action of (K1_4.P1,K2_4.P2) for the final row
    M.m[3][0] = fp2_mul(Gi[0].g00, M.m[1][0])
    M.m[3][0] = fp2_add(M.m[3][0], fp2_mul(Gi[0].g01, M.m[1][2]))

    M.m[3][1] = fp2_mul(Gi[0].g00, M.m[1][1])
    M.m[3][1] = fp2_add(M.m[3][1], fp2_mul(Gi[0].g01, M.m[1][3]))

    M.m[3][2] = fp2_mul(Gi[0].g10, M.m[1][0])
    M.m[3][2] = fp2_add(M.m[3][2], fp2_mul(Gi[0].g11, M.m[1][2]))

    M.m[3][3] = fp2_mul(Gi[0].g10, M.m[1][1])
    M.m[3][3] = fp2_add(M.m[3][3], fp2_mul(Gi[0].g11, M.m[1][3]))


    return M, True

def apply_isomorphism_general(M: basis_change_matrix_t,
                              P: ThetaPoint,
                              Pt_not_zero: bool):

    res = ThetaPoint()
    temp = ThetaPoint()

    # ---- x coordinate ----
    temp.x = fp2_mul(P.x, M.m[0][0])
    temp.x = fp2_add(temp.x, fp2_mul(P.y, M.m[0][1]))
    temp.x = fp2_add(temp.x, fp2_mul(P.z, M.m[0][2]))

    # ---- y coordinate ----
    temp.y = fp2_mul(P.x, M.m[1][0])
    temp.y = fp2_add(temp.y, fp2_mul(P.y, M.m[1][1]))
    temp.y = fp2_add(temp.y, fp2_mul(P.z, M.m[1][2]))

    # ---- z coordinate ----
    temp.z = fp2_mul(P.x, M.m[2][0])
    temp.z = fp2_add(temp.z, fp2_mul(P.y, M.m[2][1]))
    temp.z = fp2_add(temp.z, fp2_mul(P.z, M.m[2][2]))

    # ---- t coordinate ----
    temp.t = fp2_mul(P.x, M.m[3][0])
    temp.t = fp2_add(temp.t, fp2_mul(P.y, M.m[3][1]))
    temp.t = fp2_add(temp.t, fp2_mul(P.z, M.m[3][2]))

    # ---- extra terms if P.t != 0 ----
    if Pt_not_zero:
        temp.x = fp2_add(temp.x, fp2_mul(P.t, M.m[0][3]))
        temp.y = fp2_add(temp.y, fp2_mul(P.t, M.m[1][3]))
        temp.z = fp2_add(temp.z, fp2_mul(P.t, M.m[2][3]))
        temp.t = fp2_add(temp.t, fp2_mul(P.t, M.m[3][3]))

    # ---- copy to result (deep copy semantics) ----
    res.x = fp2_copy(temp.x)
    res.y = fp2_copy(temp.y)
    res.z = fp2_copy(temp.z)
    res.t = fp2_copy(temp.t)

    return res


def apply_isomorphism(M: basis_change_matrix_t, 
                      P: ThetaPoint):

    res = ThetaPoint()
    res = apply_isomorphism_general(M, P, True)

    return res




def base_change( phi: theta_gluing_t, 
                 T: theta_couple_point_t):

    out = ThetaPoint()
    null_point = ThetaPoint()

    # null_point = (a : b : c : d)
    # a = P1.x P2.x, b = P1.x P2.z, c = P1.z P2.x, d = P1.z P2.z
    null_point.x = fp2_mul(T.P1.x, T.P2.x)
    null_point.y = fp2_mul(T.P1.x, T.P2.z)
    null_point.z = fp2_mul(T.P2.x, T.P1.z)
    null_point.t = fp2_mul(T.P1.z, T.P2.z)

    #Apply the basis change
    out = apply_isomorphism(phi.M, null_point)

    return out


def gluing_compute( E12: theta_couple_curve_t,
                    xyK1_8: theta_couple_jac_point_t,
                    xyK2_8: theta_couple_jac_point_t,
                    verify: bool):

    out = theta_gluing_t()
    out.xyK1_8 = xyK1_8
    out.domain = E12

    #Given points in E[8] x E[8] we need the four torsion below
    xyK1_4 = theta_couple_jac_point_t()
    xyK2_4 = theta_couple_jac_point_t()

    xyK1_4 = double_couple_jac_point(xyK1_8, E12)
    xyK2_4 = double_couple_jac_point(xyK2_8, E12)

    #Convert from (X:Y:Z) coordinates to (X:Z)
    K1_8 = theta_couple_point_t() 
    K2_8 = theta_couple_point_t()
    K1_4 = theta_couple_point_t() 
    K2_4 = theta_couple_point_t()

    K1_8 = couple_jac_to_xz(xyK1_8)
    K2_8 = couple_jac_to_xz(xyK2_8)
    K1_4 = couple_jac_to_xz(xyK1_4)
    K2_4 = couple_jac_to_xz(xyK2_4)

    # Set the basis change matrix, if we have not been given a valid K[8] for this computation
    # gluing_change_of_basis will detect this and return 0
    out.M, ret = gluing_change_of_basis(K1_4, K2_4, E12)

    if not ret:
        print("gluing failed as kernel does not have correct order")
        return 0

    # apply the base change to the kernel
    TT1 = ThetaPoint()
    TT2 = ThetaPoint()

    TT1 = base_change(out, K1_8)
    TT2 = base_change(out, K2_8)

    # compute the codomain
    TT1 = to_squared_theta(TT1)
    TT2 = to_squared_theta(TT2)

    # If the kernel is well formed then TT1.t and TT2.t are zero
    # if they are not, we exit early as the signature we are validating
    # is probably malformed
    if (not (fp2_is_zero(TT1.t) and fp2_is_zero(TT2.t))):
        print("gluing failed TT1.t or TT2.t is not zero")
        return 0
    
    # Test our projective factors are non zero
    if (fp2_is_zero(TT1.x) or fp2_is_zero(TT2.x) or fp2_is_zero(TT1.y) or fp2_is_zero(TT2.z) or fp2_is_zero(TT1.z)):
        return 0; # invalid input

    # Projective factor: Ax
    out.codomain.x = fp2_mul(TT1.x, TT2.x)
    out.codomain.y = fp2_mul(TT1.y, TT2.x)
    out.codomain.z = fp2_mul(TT1.x, TT2.z)
    out.codomain.t = fp2_set_zero()
    # Projective factor: ABCxz

    ret_val_x = out.codomain.x
    ret_val_y = out.codomain.y
    ret_val_z = out.codomain.z
    ret_val_t = out.codomain.t

    out.precomputation = ThetaPoint()

    out.precomputation.x = fp2_mul(TT1.y, TT2.z)
    out.precomputation.y = fp2_copy(ret_val_z)
    out.precomputation.z = fp2_copy(ret_val_y)
    out.precomputation.t = fp2_set_zero()

    # Compute the two components of phi(K1_8) = (x:x:y:y).
    out.imageK1_8.x = fp2_mul(TT1.x, out.precomputation.x)
    out.imageK1_8.y = fp2_mul(TT1.z, out.precomputation.z)

    # If K1_8 and K2_8 are our 8-torsion points, this ensures that the
    # 4-torsion points [2]K1_8 and [2]K2_8 are isotropic.
    if (verify):

        t1 = fp2_mul(TT1.y, out.precomputation.y)
        if ( not fp2_is_equal(out.imageK1_8.x, t1)):
            return 0
        

        t1 = fp2_mul(TT2.x, out.precomputation.x)
        t2 = fp2_mul(TT2.z, out.precomputation.z)
        if (not fp2_is_equal(t2, t1)):
            return 0

    # compute the final codomain
    out.codomain = hadamard_sqisign(out.codomain)
    return out, True



def gluing_eval_point_special_case( P: theta_couple_point_t, phi: theta_gluing_t):
    
    image   = ThetaPoint()
    T       = ThetaPoint()

    # Apply the basis change
    T = base_change(phi, P)

    # Apply the to_squared_theta transform
    T = to_squared_theta(T)

    # This coordinate should always be 0 in a gluing because D=0.
    # If this is not the case, something went very wrong, so reject
    if not fp2_is_zero(T.t):
        return False

    # Compute (x, y, z, t)
    image.x = fp2_mul(T.x, phi.precomputation.x)
    image.y = fp2_mul(T.y, phi.precomputation.y)
    image.z = fp2_mul(T.z, phi.precomputation.z)
    image.t = fp2_set_zero()

    image = hadamard_sqisign(image)
    return image, True

def jac_to_xz_add_components(P: jac_point_t, 
                             Q: jac_point_t, 
                             AC: ECCurve)->add_components_t:

    # Take P and Q in E distinct, two jac_point_t, return three components u,v and w in Fp2 such
    # that the xz coordinates of P+Q are (u-v:w) and of P-Q are (u+v:w)

    #fp2_t t0, t1, t2, t3, t4, t5, t6;

    add_comp = add_components_t()

    t0          = fp2_sqr(P.z);             # t0 = z1^2
    t1          = fp2_sqr(Q.z);             # t1 = z2^2
    t2          = fp2_mul(P.x, t1);        # t2 = x1z2^2
    t3          = fp2_mul(t0, Q.x);        # t3 = z1^2x2
    t4          = fp2_mul(P.y, Q.z);      # t4 = y1z2
    t4          = fp2_mul(t4, t1);          # t4 = y1z2^3
    t5          = fp2_mul(P.z, Q.y);      # t5 = z1y2
    t5          = fp2_mul(t5, t0);          # t5 = z1^3y2
    t0          = fp2_mul(t0, t1);          # t0 = (z1z2)^2
    t6          = fp2_mul(t4, t5);          # t6 = (z1z_2)^3y1y2
    add_comp.v  = fp2_add(t6, t6); # v  = 2(z1z_2)^3y1y2
    t4          = fp2_sqr(t4);               # t4 = y1^2z2^6
    t5          = fp2_sqr(t5);               # t5 = z1^6y_2^2
    t4          = fp2_add(t4, t5);          # t4 = z1^6y_2^2 + y1^2z2^6
    t5          = fp2_add(t2, t3);          # t5 = x1z2^2 +z_1^2x2
    t6          = fp2_add(t3, t3);          # t6 = 2z_1^2x2
    t6          = fp2_sub(t5, t6);          # t6 = lambda = x1z2^2 - z_1^2x2
    t6          = fp2_sqr(t6);               # t6 = lambda^2 = (x1z2^2 - z_1^2x2)^2
    t1          = fp2_mul(AC.A, t0);       # t1 = A*(z1z2)^2
    t1          = fp2_add(t5, t1);          # t1 = gamma =A*(z1z2)^2 + x1z2^2 +z_1^2x2
    t1          = fp2_mul(t1, t6);          # t1 = gamma*lambda^2
    add_comp.u  = fp2_sub(t4, t1); # u  = z1^6y_2^2 + y1^2z2^6 - gamma*lambda^2
    add_comp.w  = fp2_mul(t6, t0); # w  = (z1z2)^2(lambda)^2

    return add_comp


def gluing_eval_point( P: theta_couple_jac_point_t, 
                       phi: theta_gluing_t):

    image = ThetaPoint()
    T1 = ThetaPoint()
    T2 = ThetaPoint()
    add_comp1 = add_components_t()
    add_comp2 = add_components_t()

    # Compute the cross addition components of P1+Q1 and P2+Q2
    add_comp1 = jac_to_xz_add_components(P.P1, phi.xyK1_8.P1, phi.domain.E1)
    add_comp2 = jac_to_xz_add_components(P.P2, phi.xyK1_8.P2, phi.domain.E2)

    # Compute T1 and T2 derived from the cross addition components.
    T1.x = fp2_mul(add_comp1.u, add_comp2.u);   # T1x = u1u2
    T2.t = fp2_mul(add_comp1.v, add_comp2.v);   # T2t = v1v2
    T1.x = fp2_add(T1.x, T2.t);                 # T1x = u1u2 + v1v2
    T1.y = fp2_mul(add_comp1.u, add_comp2.w);   # T1y = u1w2
    T1.z = fp2_mul(add_comp1.w, add_comp2.u);   # T1z = w1u2
    T1.t = fp2_mul(add_comp1.w, add_comp2.w);   # T1t = w1w2
    T2.x = fp2_add(add_comp1.u, add_comp1.v);   # T2x = (u1+v1)
    T2.y = fp2_add(add_comp2.u, add_comp2.v);   # T2y = (u2+v2)
    T2.x = fp2_mul(T2.x, T2.y);                 # T2x = (u1+v1)(u2+v2)
    T2.x = fp2_sub(T2.x, T1.x);                 # T1x = v1u2 + u1v2
    T2.y = fp2_mul(add_comp1.v, add_comp2.w);   # T2y = v1w2
    T2.z = fp2_mul(add_comp1.w, add_comp2.v);   # T2z = w1v2
    T2.t = fp2_set_zero();                      # T2t = 0

    # Apply the basis change and compute their respective square
    # theta(P+Q) = M.T1 - M.T2 and theta(P-Q) = M.T1 + M.T2
    T1 = apply_isomorphism_general(phi.M, T1, True)
    T2 = apply_isomorphism_general(phi.M, T2, False)
    T1 = pointwise_square(T1)
    T2 = pointwise_square(T2)

    # the difference between the two is therefore theta(P+Q)theta(P-Q)
    # whose hadamard transform is then the product of the dual
    # theta_points of phi(P) and phi(Q).
    T1.x = fp2_sub(T1.x, T2.x)
    T1.y = fp2_sub(T1.y, T2.y)
    T1.z = fp2_sub(T1.z, T2.z)
    T1.t = fp2_sub(T1.t, T2.t)
    T1 = hadamard_sqisign(T1)

    # Compute (x, y, z, t)
    # As imageK1_8 = (x:x:y:y), its inverse is (y:y:x:x).
    image.x = fp2_mul(T1.x, phi.imageK1_8.y)
    image.y = fp2_mul(T1.y, phi.imageK1_8.y)
    image.z = fp2_mul(T1.z, phi.imageK1_8.x)
    image.t = fp2_mul(T1.t, phi.imageK1_8.x)

    image = hadamard_sqisign(image)

    return image




def gluing_eval_basis( xyT1: theta_couple_jac_point_t,
                       xyT2: theta_couple_jac_point_t,
                       phi: theta_gluing_t):

    image1 = ThetaPoint()
    image2 = ThetaPoint()
    image1 = gluing_eval_point(xyT1, phi)
    image2 = gluing_eval_point(xyT2, phi)

    return image1, image2



def theta_isogeny_compute(
    A: ThetaStructure,
    T1_8: ThetaPoint,
    T2_8: ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
    verify: bool) -> ThetaIsogeny:

    out = ThetaIsogeny()
    # Assign top-level isogeny attributes
    out.hadamard_bool_1 = hadamard_bool_1
    out.hadamard_bool_2 = hadamard_bool_2
    out.domain = A
    out.T1_8 = T1_8
    out.T2_8 = T2_8
    out.codomain.precomputation = False

    # Local theta points (temporary variables)
    TT1 = ThetaPoint(
        x=Fp2(0, 0),
        y=Fp2(0, 0),
        z=Fp2(0, 0),
        t=Fp2(0, 0)
    )

    TT2 = ThetaPoint(
        x=Fp2(0, 0),
        y=Fp2(0, 0),
        z=Fp2(0, 0),
        t=Fp2(0, 0)
    )

    # If Hadamard transformation enabled
    if hadamard_bool_1:
        TT1 = hadamard_sqisign(T1_8) # 4 add + 4 sub
        TT1 = to_squared_theta(TT1)  # 2 add + 1 sub + 2 mul

        TT2 = hadamard_sqisign(T2_8) # 4 add + 4 sub
        TT2 = to_squared_theta(TT2) # 2 add + 1 sub + 2 mul
    else:
        TT1 = to_squared_theta(T1_8) # 2 add + 1 sub + 2 mul
        TT2 = to_squared_theta(T2_8) # 2 add + 1 sub + 2 mul

    if (
        fp2_is_zero(TT2.x)
        or fp2_is_zero(TT2.y)
        or fp2_is_zero(TT2.z)
        or fp2_is_zero(TT2.t)
        or fp2_is_zero(TT1.x)
        or fp2_is_zero(TT1.y)
    ):
        return False  # corresponds to 'return 0;' in C
    

    t1 = fp2_mul(TT1.x, TT2.y) # 2 add 3 sub 3 mul
    t2 = fp2_mul(TT1.y, TT2.x) # 2 add 3 sub 3 mul
    out.codomain.null_point.x = fp2_mul(TT2.x, t1)
    out.codomain.null_point.y = fp2_mul(TT2.y, t2)
    out.codomain.null_point.z = fp2_mul(TT2.z, t1)
    out.codomain.null_point.t = fp2_mul(TT2.t, t2)

    t3 = fp2_mul(TT2.z, TT2.t) # 2 add 3 sub 3 mul
    out.precomputation.x = fp2_mul(t3, TT1.y)
    out.precomputation.y = fp2_mul(t3, TT1.x)
    out.precomputation.z = fp2_copy(out.codomain.null_point.t)
    out.precomputation.t = fp2_copy(out.codomain.null_point.z)
    

    if verify:
        # (1) TT1.x * out.precomputation.x == TT1.y * out.precomputation.y
        t1 = fp2_mul(TT1.x, out.precomputation.x)
        t2 = fp2_mul(TT1.y, out.precomputation.y)
        if not fp2_is_equal(t1, t2):
            return False
        
        # (2) TT1.z * out.precomputation.z == TT1.t * out.precomputation.t
        t1 = fp2_mul(TT1.z, out.precomputation.z)
        t2 = fp2_mul(TT1.t, out.precomputation.t)
        if not fp2_is_equal(t1, t2):
            return False

        # (3) TT2.x * out.precomputation.x == TT2.z * out.precomputation.z
        t1 = fp2_mul(TT2.x, out.precomputation.x)
        t2 = fp2_mul(TT2.z, out.precomputation.z)
        if not fp2_is_equal(t1, t2):
            return False

        # (4) TT2.y * out.precomputation.y == TT2.t * out.precomputation.t
        t1 = fp2_mul(TT2.y, out.precomputation.y)
        t2 = fp2_mul(TT2.t, out.precomputation.t)
        if not fp2_is_equal(t1, t2):
            return False

    if hadamard_bool_2:
        out.codomain.null_point = hadamard_sqisign(out.codomain.null_point) # 4 add + 4 sub

    return out, True

def double_point(A: ThetaStructure, inp: ThetaPoint) -> ThetaPoint:
    """
    Python equivalent of:
        void double_point(theta_point_t *out, theta_structure_t *A, const theta_point_t *in)
    Performs a single theta-point doubling with precomputation and fp2 operations.
    """
    # --- Step 1: Move to squared theta coordinates
    out = to_squared_theta(inp)

    # --- Step 2: Square each coordinate
    out.x = fp2_sqr(out.x)
    out.y = fp2_sqr(out.y)
    out.z = fp2_sqr(out.z)
    out.t = fp2_sqr(out.t)

    # --- Step 3: Ensure precomputation
    if not A.precomputation:
        A = theta_precomputation(A)

    # --- Step 4: Multiply with A’s precomputed constants
    out.x = fp2_mul(out.x, A.YZT0)
    out.y = fp2_mul(out.y, A.XZT0)
    out.z = fp2_mul(out.z, A.XYT0)
    out.t = fp2_mul(out.t, A.XYZ0)

    # --- Step 5: Hadamard transform
    out = hadamard_sqisign(out)

    # --- Step 6: Final multiplication with lowercase constants
    out.x = fp2_mul(out.x, A.yzt0)
    out.y = fp2_mul(out.y, A.xzt0)
    out.z = fp2_mul(out.z, A.xyt0)
    out.t = fp2_mul(out.t, A.xyz0)

    return out, A


def double_iter(A: ThetaStructure, inp: ThetaPoint, exp: int) -> ThetaPoint:
    """
    Python equivalent of:
        void double_iter(theta_point_t *out, theta_structure_t *A, const theta_point_t *in, int exp)
    Performs repeated theta-doubling.
    """
    A_new = ThetaStructure()
    A_new = A
    # Base case — exp == 0
    if exp == 0:
        return ThetaPoint(
            x=Fp2(inp.x.re, inp.x.im),
            y=Fp2(inp.y.re, inp.y.im),
            z=Fp2(inp.z.re, inp.z.im),
            t=Fp2(inp.t.re, inp.t.im)
        )

    # Perform the first doubling
    out, A_new = double_point(A_new, inp)

    # Repeat exp−1 additional doublings
    for _ in range(1, exp):
        out, A_new = double_point(A_new, out)

    return out, A_new





def theta_isogeny_eval(phi: ThetaIsogeny, P: ThetaPoint):
    """
    Python equivalent of:
        void theta_isogeny_eval(theta_point_t *out, const theta_isogeny_t *phi, const theta_point_t *P)
    """

    out = ThetaPoint()

    # Step 1: Apply hadamard or direct square transformation
    if phi.hadamard_bool_1:
        calc1 = hadamard_sqisign(P)
        calc1 = to_squared_theta(calc1)
    else:
        calc1 = to_squared_theta(P)

    # Step 2: Coordinate-wise fp2 multiplications with precomputation
    out.x = fp2_mul(calc1.x, phi.precomputation.x)
    out.y = fp2_mul(calc1.y, phi.precomputation.y)
    out.z = fp2_mul(calc1.z, phi.precomputation.z)
    out.t = fp2_mul(calc1.t, phi.precomputation.t)

    # Step 3: Optional hadamard transformation (based on hadamard_bool_2)
    if phi.hadamard_bool_2:
        calc2 = hadamard_sqisign(out)
        out = calc2

    return out

def choose_index_theta_point(ind: int, T: ThetaPoint):

    res = Fp2()
    src = Fp2()
    r = ind % 4

    if r == 0:
        src = T.x
    elif r == 1:
        src = T.y
    elif r == 2:
        src = T.z
    elif r == 3:
        src = T.t
    else:
        raise AssertionError("unreachable")


    res = fp2_copy(src)

    return res

def select_base_change_matrix( M1: basis_change_matrix_t, M2: basis_change_matrix_t, option: int)->basis_change_matrix_t:
    """
    If option == 0:
        M <- M1
    else:
        M <- FP2_CONSTANTS[M2]
    """
    
    M = basis_change_matrix_init_zero()
    
    for i in range(4):
        for j in range(4):
            M.m[i][j] = fp2_select(
                M1.m[i][j],
                FP2_CONSTANTS[M2[i][j]],
                option
            )

    return M

def splitting_compute(A: ThetaStructure,
                      zero_index: int,
                      randomize: bool) -> bool:
    
    out = theta_splitting_t()
    # init
    count = 0
    U_cst = fp2_set_zero()

    # zero the splitting matrix
    out.M = basis_change_matrix_init_zero()

    # enumerate through all indices
    for i in range(10):
        U_cst = fp2_set_zero()

        for t in range(4):
            # iterate through the null point
            t2 = choose_index_theta_point(t, A.null_point)
            t1 = choose_index_theta_point(t ^ EVEN_INDEX[i][1], A.null_point)

            # t1 = t1 * t2
            t1 = fp2_mul(t1, t2)

            # ctl logic derived from CHI_EVAL
            ctl = CHI_EVAL[EVEN_INDEX[i][0]][t] >> 1
            assert ctl in (0, -1, 0xFFFFFFFF)

            # conditional sign
            t2 = fp2_neg(t1)
            t1 = fp2_select(t1, t2, ctl)

            # accumulate
            U_cst = fp2_add(U_cst, t1)

        # if U_cst == 0, update splitting matrix
        ctl = fp2_is_zero(U_cst)
        if ctl == False:
            ctl = 0
        else:
            ctl = 1
        count += ctl
        out.M = select_base_change_matrix(out.M, SPLITTING_TRANSFORMS[i], ctl)

        # extra consistency check
        if zero_index != -1 and i == zero_index and not ctl:
            return False


    # apply isomorphism to ensure compatibility
    out.B = ThetaStructure()
    out.B.null_point = apply_isomorphism(
        out.M,
        A.null_point
    )

    # success iff exactly one zero was found
    return count == 1, out


def is_product_theta_point(P: ThetaPoint):    
    t1 = fp2_mul( P.x, P.t)
    t2 = fp2_mul( P.y, P.z)
    return fp2_is_equal(t1, t2)


def theta_product_structure_to_elliptic_product(A: ThetaStructure):
    # fp2_t xx, yy;

    E12 = theta_couple_curve_t(E1=ec_curve_init(), E2=ec_curve_init)
    # This should be true from our computations in splitting_compute
    # but still check this for sanity
    ret = is_product_theta_point(A.null_point)
    if not ret:
        return False
    
    E12.E1 = ec_curve_init()
    E12.E2 = ec_curve_init()

    # A valid elliptic theta null point has no zero coordinate
    if (fp2_is_zero(A.null_point.x) or fp2_is_zero(A.null_point.y) or fp2_is_zero(A.null_point.z)):
        return 0

    # xx = x², yy = y²
    xx = fp2_sqr(A.null_point.x)
    yy = fp2_sqr(A.null_point.y)
    # xx = x^4, yy = y^4
    xx = fp2_sqr(xx)
    yy = fp2_sqr(yy)

    # A2 = -2(x^4+y^4)/(x^4-y^4)
    E12.E2.A = fp2_add(xx, yy)
    E12.E2.C = fp2_sub(xx, yy)
    E12.E2.A = fp2_add(E12.E2.A, E12.E2.A)
    E12.E2.A = fp2_neg(E12.E2.A)

    # same with x,z
    xx = fp2_sqr(A.null_point.x)
    yy = fp2_sqr(A.null_point.z)
    xx = fp2_sqr(xx)
    yy = fp2_sqr(yy)

    # A1 = -2(x^4+z^4)/(x^4-z^4)
    E12.E1.A = fp2_add(xx, yy)
    E12.E1.C = fp2_sub(xx, yy)
    E12.E1.A = fp2_add(E12.E1.A, E12.E1.A)
    E12.E1.A = fp2_neg(E12.E1.A)

    if (fp2_is_zero(E12.E1.C) or fp2_is_zero(E12.E2.C)):
        return False

    return E12,True


def theta_point_to_montgomery_point(P12_in: theta_couple_point_t, P: ThetaPoint, A: ThetaStructure):
    P12 = theta_couple_point_t(P1=ECPoint(x=P12_in.P1.x, z=P12_in.P1.z), P2=ECPoint(x=P12_in.P2.x, z=P12_in.P2.z))

    ret = is_product_theta_point(P)
    if not ret:
        return False

    x = fp2_copy(P.x)
    z = fp2_copy(P.y)
    if (fp2_is_zero(x) and fp2_is_zero(z)):
        x = fp2_copy(P.z)
        z = fp2_copy(P.t)
    
    if (fp2_is_zero(x) and fp2_is_zero(z)):
        return 0 # at this point P=(0:0:0:0) so is invalid
    
    # P2.X = A.null_point.y * P.x + A.null_point.x * P.y
    # P2.Z = - A.null_point.y * P.x + A.null_point.x * P.y
    P12.P2.x = fp2_mul(A.null_point.y, x)
    temp     = fp2_mul(A.null_point.x, z)
    P12.P2.z = fp2_sub(temp, P12.P2.x)
    P12.P2.x = fp2_add(P12.P2.x, temp)

    x = fp2_copy(P.x)
    z = fp2_copy(P.z)
    if (fp2_is_zero(x) and fp2_is_zero(z)):
        x = fp2_copy(P.y)
        z = fp2_copy(P.t)
    
    # P1.X = A.null_point.z * P.x + A.null_point.x * P.z
    # P1.Z = -A.null_point.z * P.x + A.null_point.x * P.z
    P12.P1.x    = fp2_mul(A.null_point.z, x)
    temp        = fp2_mul(A.null_point.x, z)
    P12.P1.z   = fp2_sub(temp, P12.P1.x)
    P12.P1.x   = fp2_add(P12.P1.x, temp)
    return P12 ,True


def _theta_chain_compute_impl( n: int,
                               E12_in: theta_couple_curve_t,
                               ker: theta_kernel_couple_points_t,
                               extra_torsion: bool,
                               numP: int,
                               verify: bool,
                               randomize: bool):

    E12 = theta_couple_curve_t(E1=E12_in.E1, E2=E12_in.E2)
    E34 = theta_couple_curve_t(E1=ec_curve_init, E2=ec_curve_init)
    #P12 = list([theta_couple_point_t for i in range()])
    theta = ThetaStructure()

    #lift the basis
    xyT1 = theta_couple_jac_point_t(P1=jac_point_t, P2=jac_point_t)
    xyT2 = theta_couple_jac_point_t(P1=jac_point_t, P2=jac_point_t)

    bas1 = ECBasis(
        P=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        Q=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        PmQ=ECPoint(Fp2(0, 0), Fp2(0, 0)),
    )
    bas1.P = ker.T1.P1
    bas1.Q = ker.T2.P1
    bas1.PmQ = ker.T1m2.P1
    bas2 =  ECBasis(
        P=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        Q=ECPoint(Fp2(0, 0), Fp2(0, 0)),
        PmQ=ECPoint(Fp2(0, 0), Fp2(0, 0)),
    )
    bas2.P = ker.T1.P2
    bas2.Q = ker.T2.P2
    bas2.PmQ = ker.T1m2.P2
    xyT1.P1, xyT2.P1, ret = lift_basis(bas1, E12.E1)
    if not ret:
        return 0
    xyT1.P2, xyT2.P2, ret2 = lift_basis(bas2, E12.E2)
    if not ret2:
        return 0
    extra = HD_extra_torsion * extra_torsion
    if numP != None:
        size_nump = numP
    else:
        size_nump = 1
    

    pts = [ThetaPoint() for i in range(size_nump)]

    space = 1
    #for i = 1; i < n; i *= 2
    i = 1
    while i < n:
        space = space + 1
        i <<= 1   # same as i *= 2

    todo = [0 for i in range(space)]
    todo[0] = n - 2 + extra

    current = 0

    #kernel points for the gluing isogeny
    jacQ1 = [theta_couple_jac_point_t(P1=jac_point_t(x=Fp2(0, 0), y=Fp2(0, 0), z=Fp2(0, 0)), P2=jac_point_t(x=Fp2(0, 0), y=Fp2(0, 0), z=Fp2(0, 0))) for i in range(space)]
    jacQ2 = [theta_couple_jac_point_t(P1=jac_point_t(x=Fp2(0, 0), y=Fp2(0, 0), z=Fp2(0, 0)), P2=jac_point_t(x=Fp2(0, 0), y=Fp2(0, 0), z=Fp2(0, 0))) for i in range(space)]
    jacQ1[0] = xyT1
    jacQ2[0] = xyT2
    while todo[current] != 1:
        assert todo[current] >= 2 
        current = current + 1
        assert current < space 
        # the gluing isogeny is quite a bit more expensive than the others,
        # so we adjust the usual splitting rule here a little bit: towards
        # the end of the doubling chain it will be cheaper to recompute the
        # doublings after evaluation than to push the intermediate points.
        if todo[current - 1] >= 16:
            num_dbls =  todo[current - 1] // 2 
        else:
            num_dbls = todo[current - 1] - 1
        #num_dbls = todo[current - 1] >= 16 ? todo[current - 1] / 2 : todo[current - 1] - 1;
        assert num_dbls > 0 and num_dbls < todo[current - 1]
        jacQ1[current] = double_couple_jac_point_iter(num_dbls, jacQ1[current - 1], E12)
        jacQ2[current] = double_couple_jac_point_iter(num_dbls, jacQ2[current - 1], E12)
        todo[current] = todo[current - 1] - num_dbls
    

    # kernel points for the remaining isogeny steps
    thetaQ1 = [ThetaPoint() for i in range(space)]
    thetaQ2 = [ThetaPoint() for i in range(space)]

    #the gluing step
    first_step = theta_gluing_t()
    
    assert todo[current] == 1

    # compute the gluing isogeny
    first_step, ret = gluing_compute(E12, jacQ1[current], jacQ2[current], verify)
    if not ret:
        return False

    # evaluate
    for j in range(numP):
        assert ec_is_zero(P12[j].P1) or ec_is_zero(P12[j].P2)
        pts[j], ret = gluing_eval_point_special_case( P12[j], first_step)
        if not ret:
            return 0
    

    #push kernel points through gluing isogeny
    for j in range(current):
        thetaQ1[j], thetaQ2[j] = gluing_eval_basis(jacQ1[j], jacQ2[j], first_step)
        todo[j] = todo[j] - 1
    


    current = current - 1
    

    # set-up the theta_structure for the first codomain
    theta.null_point = first_step.codomain
    theta.precomputation = 0
    theta = theta_precomputation(theta)

    step = theta_isogeny_t()

    # and now we do the remaining steps
    i = 1
    while current >= 0 and todo[current] != 0:
        assert current < space 
        while todo[current] != 1:
            assert todo[current] >= 2
            current = current + 1
            assert current < space
            num_dbls = todo[current - 1] // 2
            assert num_dbls and num_dbls < todo[current - 1]
            thetaQ1[current], theta = double_iter(theta, thetaQ1[current - 1], num_dbls)
            thetaQ2[current], theta = double_iter(theta, thetaQ2[current - 1], num_dbls)
            todo[current] = todo[current - 1] - num_dbls


        #computing the next step

        step = ThetaIsogeny()
        
        if i == n - 2: # penultimate step
            step, ret = theta_isogeny_compute(theta, thetaQ1[current], thetaQ2[current], 0, 0, verify)
        elif i == n - 1: # ultimate step
            step, ret = theta_isogeny_compute(theta, thetaQ1[current], thetaQ2[current], 1, 0, False)
        else:
            step, ret = theta_isogeny_compute(theta, thetaQ1[current], thetaQ2[current], 0, 1, verify)
        

        if not ret:
            return 0

        for j in range(numP):
            pts[j] = theta_isogeny_eval(step, pts[j])

        # updating the codomain
        theta = step.codomain

        # pushing the kernel
        assert todo[current] == 1
        for j in range(current):
            thetaQ1[j] = theta_isogeny_eval(step, thetaQ1[j])
            thetaQ2[j] = theta_isogeny_eval(step, thetaQ2[j])
            assert todo[j] > 0
            todo[j] = todo[j] - 1
        current = current - 1
        i = i + 1
    assert current == -1

    if not extra_torsion:
        if n >= 3:
            # in the last step we've skipped pushing the kernel since current was == 0, let's do it now
            thetaQ1[0] = theta_isogeny_eval(step, thetaQ1[0])
            thetaQ2[0] = theta_isogeny_eval(step, thetaQ2[0])

        # penultimate step
        step = theta_isogeny_compute_4(theta, thetaQ1[0], thetaQ2[0], 0, 0)
        for j in range(numP):
            pts[j] = theta_isogeny_eval(step, pts[j])
        theta = step.codomain
        thetaQ1[0] = theta_isogeny_eval( step, thetaQ1[0])
        thetaQ2[0] = theta_isogeny_eval( step, thetaQ2[0])

        # ultimate step
        step = theta_isogeny_compute_2(theta, thetaQ1[0], thetaQ2[0], 1, 0)
        for j in range(numP):
            pts[j] = theta_isogeny_eval( step, pts[j])
        theta = step.codomain
    

    # final splitting step
    last_step = theta_splitting_t()

    if extra_torsion:
        split_in1 = 8
    else:
        split_in1 = -1

    is_split, last_step = splitting_compute(theta, split_in1 , randomize)

    if not is_split:
        print("kernel did not generate an isogeny between elliptic products")
        return 0
    
    E34, ret = theta_product_structure_to_elliptic_product(last_step.B)
    if not ret: 
        return 0

    # evaluate
    for j in range(numP):
        pts[j] = apply_isomorphism(last_step.M, pts[j])
        P12[j], ret = theta_point_to_montgomery_point( pts[j], last_step.B)
        if not ret:
            return False
    

    return True, E12 , E34



def theta_chain_compute_and_eval_verify(n: int,
                                    E12: theta_couple_curve_t,
                                    ker:theta_kernel_couple_points_t,
                                    extra_torsion: bool,
                                    numP: int):

    return _theta_chain_compute_impl(n, E12, ker, extra_torsion, numP, True, False)



