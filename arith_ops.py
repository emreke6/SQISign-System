from dataclasses import dataclass, field
from typing import List, Tuple, Optional


q   = pow(2,248)*5 - 1 # 64-bit: prime q
R   = pow(2,256) # 2^256
Rv  = pow(R, -1, q) # 64-bit: R^-1 (mod q)
mu  = (-pow(q,-1,R)) % R # 64-bit: (−q)^−1 (mod R)

ONE = 0x0100000000000000000000000000000000000000000000000000000000000033

SECURITY_BITS = 128
HASH_ITERATIONS = 64
FP2_ENCODED_BYTES = 64

RADIX = 64
LOG2RADIX = 6




R_INV = pow(R, -1, q)   # modular inverse

# ---------- Base Fp2 ----------
@dataclass
class Fp2:
    re: int = 0
    im: int = 0

@dataclass
class ECPoint:
    x: Fp2
    z: Fp2

@dataclass
class ECCurve:
    A: Fp2
    C: Fp2
    A24: Optional[ECPoint]
    is_A24_computed_and_normalized: bool


# ---------- Theta Point ----------
@dataclass
class ThetaPoint:
    x: Fp2 = field(default_factory=Fp2)
    y: Fp2 = field(default_factory=Fp2)
    z: Fp2 = field(default_factory=Fp2)
    t: Fp2 = field(default_factory=Fp2)


# ---------- Theta Structure ----------
@dataclass
class ThetaStructure:
    null_point: ThetaPoint = field(default_factory=ThetaPoint)
    precomputation: bool = False

    # Eight precomputed fp2 values
    XYZ0: Fp2 = field(default_factory=Fp2)
    YZT0: Fp2 = field(default_factory=Fp2)
    XZT0: Fp2 = field(default_factory=Fp2)
    XYT0: Fp2 = field(default_factory=Fp2)
    xyz0: Fp2 = field(default_factory=Fp2)
    yzt0: Fp2 = field(default_factory=Fp2)
    xzt0: Fp2 = field(default_factory=Fp2)
    xyt0: Fp2 = field(default_factory=Fp2)


# ---------- Theta Isogeny ----------
@dataclass
class ThetaIsogeny:
    T1_8: ThetaPoint = field(default_factory=ThetaPoint)
    T2_8: ThetaPoint = field(default_factory=ThetaPoint)
    hadamard_bool_1: bool = False
    hadamard_bool_2: bool = False
    domain: ThetaStructure = field(default_factory=ThetaStructure)
    precomputation: ThetaPoint = field(default_factory=ThetaPoint)
    codomain: ThetaStructure = field(default_factory=ThetaStructure)




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



def montgomery(C, R, q):
    m = ((C%(R)) * (mu)) % R
    u = C + m*(q)
    u = u >> 256 # u/R
    D = u
    return D


def interleaved_montgomery(C, q):
    """
    64-bit interleaved Montgomery reduction for 256-bit modulus q.
    C: integer (512-bit product)
    q: integer (modulus)
    returns: C * R^-1 mod q, where R = 2^256
    """
    R_word = 1 << 64
    mask = R_word - 1

    # precompute q' = -q^{-1} mod 2^64
    q_dash = (-pow(q, -1, R_word)) & mask

    # work on 320-bit intermediate value
    T = C
    for i in range(4):  # 4 limbs of 64 bits
        m = (((T) & mask) * (q_dash)) & mask
        T = (T + m * (q)) >> 64

    print("check: ", T > q, T > 2*q)
    # conditional subtraction

    return T


def montgomery_sqisign(C, R, q):
    C_new = C % R
    C_times_mu = (C_new << 250) + (C_new << 248) + C_new # C*q
    m = C_times_mu % (pow(2,256))
    m_times_q_calc = (((m<<2) + m) << 248) - m # m*q^(-1) mod R
    u = C + m_times_q_calc
    u = u >> 256 # u/R
    # if u > q:
    #     D = u - q
    # else:
    #     D = u
    D = u
    return D


def montgomery_true(C, R, q):
    C_new = C % R
    C_times_mu = (C_new << 250) + (C_new << 248) + C_new # C*q
    m = C_times_mu % (pow(2,256))
    m_times_q_calc = (((m<<2) + m) << 248) - m # m*q^(-1) mod R
    u = C + m_times_q_calc
    u = u >> 256 # u/R
    if u > q:
        D = u - q
    else:
        D = u
    #D = u
    return D

def montgomery_sqisign_fp_sqr(C, R, q):
    C_new = C
    C_times_mu = (C_new << 250) + (C_new << 248) + C_new # C*q
    m = C_times_mu % (pow(2,256))
    m_times_q_calc = (((m<<2) + m) << 248) - m # m*q^(-1) mod R
    u = C + m_times_q_calc
    u = u >> 256 # u/R
    D = u
    return D

def modsub_sqisign(a, b, bool1):
    c = (a - b) % q

    # c2 = c + q

    # if c < 0:
    #     r1 = c2
    # else:
    #     r1 = c

    # if not bool1:
    #     if r1 < 0:
    #         r1 = r1 + q
    #     else:
    #         r1 = r1
    # else:
    #     r1 = a - b + 2 * q
    r1 = c
    return r1

MASK64  = (1 << 64) - 1
MASK256 = (1 << 256) - 1

def modadd_sqisign_exact(a, b, p):
    # force 256-bit arithmetic
    c = (a + b) & MASK256

    for _ in range(2):
        # extract top limb (bits 192..255)
        r11 = (c >> 192) & MASK64

        # extract overflow bits (top 5 bits)
        overflow = r11 >> 59

        # create mask: 0 or -1 (64-bit)
        mask = (-overflow) & MASK64

        # masked subtraction (limb-wise)
        c = (c - (mask | (mask << 64) | (mask << 128) |
                 ((p >> 192) & MASK64) * (mask >> 63) << 192)
            ) & MASK256

    return c

def modadd_sqisign(a, b, bool1):
    c = (a + b) % q
    # prev = a + b

    # # First reduction pass
    # # overflow = (c >> 251) & 0x1F  # Check top bits (assuming 256-bit)
    # # mask = (-overflow) & ((1 << 256) - 1)  # Create mask
    # # print("mask: ", bin(mask), len(bin(mask)[2:]))
    # # print("mask & q = ", bin((mask & q)))
    # # c = (c - (mask & q)) & ((1 << 256) - 1)

    # # # Second reduction pass (repeat)
    # # overflow = (c >> 251) & 0x1F
    # # mask = (-overflow) & ((1 << 256) - 1)
    # # c = (c - (mask & q)) & ((1 << 256) - 1)

    # c_mask = c >> (192+59)

    # if c_mask != 0:
    #     c = c - q

    # c_mask = c >> (192+59)

    # if c_mask != 0:
    #     c = c - q

    # if bool1:
    #     c = prev
    # else:
    #     c = c

    return c


def fp_add_sqr(a, b):
    return modadd_sqisign(a, b, True)

def fp_add(a, b):
    return modadd_sqisign(a, b, False)

def fp_sub(a, b):
    return modsub_sqisign(a, b, False)

def fp_sub_sqr(a, b):
    return modsub_sqisign(a, b, True)

def fp_is_zero(a: int, q: int) -> bool:
    """
    Check if a field element 'a' in F_p is zero.
    In Montgomery or standard representation, 'a' is considered zero if:
        a % q == 0
    """
    return (a % q) == 0

def fp2_set_zero() -> Fp2:
    return Fp2(0, 0)

def fp2_set_one() -> Fp2:
    return Fp2(re=ONE, im=0)

def fp_set_one() -> int:
    return ONE

def fp_half(a: int) -> int:
    """
    Functional equivalent of gf5248_half.
    Treated as integer division by 2.
    """
    if a % 2 == 1:
        a = ((a >> 1) + ((q+1)>>1)) % q
    else:
        a = (a >> 1) % q
    return a

def fp2_half(y: Fp2) -> Fp2:
    """
    Functional equivalent of gf5248_half.
    Treated as integer division by 2.
    """
    x = Fp2(re=0, im=0)
    x.re = fp_half(y.re)
    x.im = fp_half(y.im)
    return x

def fp2_add_one(y: Fp2) -> Fp2:
    """
    Functional equivalent of gf5248_half.
    Treated as integer division by 2.
    """
    x = Fp2(re=0, im=0)
    x.re = fp_add(y.re, ONE)
    x.im = fp_copy(y.im)
    return x


def fp_neg(a: int) -> int:
    # res = 2*q - a
    # if res < pow(2,251):
    #     return res
    # else:
    #     return res - q
    return (-a) % q

def fp2_neg(y: Fp2) -> Fp2:
    x = Fp2(re=0, im=0)
    x.re = fp_neg(y.re)
    x.im = fp_neg(y.im)
    return x


def fp2_is_one(a: Fp2) -> bool:
    return a.re % q == ONE and a.im == 0

def fp_set_small(a: int)->int:
    return a * pow(2,256,q) % q

def fp_mul_small(a: int, n: int) -> int:
    """
    Functional equivalent of fp_mul_small / gf5248_mul_small.
    Treated as regular integer multiplication.
    """
    return (a * n) % q

def fp2_mul_small(y, n):
    """
    Functional equivalent of fp2_mul_small.
    Returns y * n in Fp2.
    """
    return Fp2(
        re=fp_mul_small(y.re, n),
        im=fp_mul_small(y.im, n),
    )

def legendre_symbol(a: int) -> int:
    """
    Returns:
        0  if a == 0 (mod p)
        1  if a is a quadratic residue mod p
       -1  if a is a non-residue mod p
    """
    a %= q
    if a == 0:
        return 0

    # Euler's criterion
    ls = pow(a, (q - 1) // 2, q)
    return 1 if ls == 1 else -1

def fp_is_square(a: int) -> bool:
    """
    Python equivalent of fp_is_square.
    Returns True if a == 0 or a is a quadratic residue mod p.
    """
    ls = legendre_symbol(a)
    return ls >= 0

def fp2_is_square(x: Fp2) -> Fp2:
    """
    Functional equivalent of C fp2_is_square.

    Returns True if x is a square in Fp2 (using norm test).
    """

    t0 = fp_mul(x.re, x.re)
    t1 = fp_mul(x.im, x.im)
    t0 = fp2_add(t0, t1)

    return fp_is_square(t0)

def fp_div3(a_in: int, q: int, is_mont: bool = False) -> int:
    """
    Matches the C fp2_sqrt() including:
      - canonical fp_sqrt() output (even low bit)
      - special-case when a1 == 0: force delta = a0
      - final canonical sign rule on (re, im)
    Assumes q is an odd prime with q % 4 == 3 and i^2 = -1.
    """
    
    if is_mont:
        in_calc = a_in * R_INV % q
    else:
        in_calc = a_in

    res = (in_calc * pow(3, -1, q)) % q

    return (res * R) % q

def fp2_is_zero(a: Fp2) -> bool:
    """
    Check if an Fp2 element (a.re, a.im) is zero.
    Returns True if both components are 0 mod q.
    """
    return fp_is_zero(a.re, q) and fp_is_zero(a.im, q)

def fp2_copy(in1: Fp2) -> Fp2:
    return Fp2(in1.re, in1.im)

def fp_copy(in1):
    return in1

def fp2_is_equal(a: Fp2, b: Fp2) -> bool:
    """Check equality of two Fp2 elements."""
    return a.re % q == b.re % q and a.im % q == b.im % q

def fp2_add(a, b):
    fp_add_re = fp_add(a.re, b.re)
    fp_add_im = fp_add(a.im, b.im)

    return Fp2(re=fp_add_re, im=fp_add_im)

def fp2_sub(a, b):
    fp_sub_re = fp_sub(a.re, b.re)
    fp_sub_im = fp_sub(a.im, b.im)

    return Fp2(re=fp_sub_re, im=fp_sub_im)



def inner_gf5248_montgomery_reduce(a: int) -> int:
    """
    Equivalent to inner_gf5248_montgomery_reduce
    Input: a in Montgomery form
    Output: canonical field element
    """
    return (a * R_INV) % q


def gf5248_encode(a: int) -> bytes:
    """
    a: integer in Montgomery form
    returns: 32-byte little-endian encoding
    """
    x = inner_gf5248_montgomery_reduce(a)

    out = bytearray(32)
    for i in range(4):
        limb = (x >> (64 * i)) & ((1 << 64) - 1)
        out[8*i : 8*i + 8] = limb.to_bytes(8, "little")

    return bytes(out)

def fp_encode(a: int) -> bytes:
    return gf5248_encode(a)


def fp2_encode(a: "Fp2") -> bytes:
    return fp_encode(a.re) + fp_encode(a.im)



def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return x % m

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def fp2_inv(x: Fp2) -> Fp2:
    # t0 = a^2
    t0 = fp_mul(x.re, x.re)

    # t1 = b^2
    t1 = fp_mul(x.im, x.im)

    # t0 = a^2 + b^2
    t0 = fp_add(t0, t1)

    # t0 = 1 / (a^2 + b^2)
    t0 = (modinv(t0 * R_INV, q) * R) % q

    # re = a * t0
    re = fp_mul(x.re, t0)

    # im = -b * t0
    t0 = fp_neg(t0)
    im = fp_mul(x.im, t0)
    

    return Fp2(re=re, im=im)


def fp_inv(a:int)->int:
    return pow(a,-1,q) * pow(2, 256, q) % q


def fp_mul(a, b):
    C = a * b
    return montgomery_sqisign(C, R, q)

def fp_mul_fp_sqr(a, b):
    C = a * b
    return montgomery_sqisign_fp_sqr(C, R , q)

def fp2_mul(y, z):
    """
    Perform fp2 multiplication:
        x = y * z
    where each input is an Fp2 object with .re and .im attributes,
    and fp_add, fp_sub, fp_mul each return Fp elements.
    """
    # Step 1: t0 = (y.re + y.im)
    t0 = fp_add(y.re, y.im)

    # Step 2: t1 = (z.re + z.im)
    t1 = fp_add(z.re, z.im)

    # Step 3: t0 = t0 * t1
    t0 = fp_mul(t0, t1)

    # Step 4: t1 = (y.im * z.im)
    t1 = fp_mul(y.im, z.im)

    # Step 5: x.re = (y.re * z.re)
    x_re = fp_mul(y.re, z.re)

    # Step 6: x.im = t0 - t1
    x_im = fp_sub(t0, t1)

    # Step 7: x.im = x.im - x.re
    x_im = fp_sub(x_im, x_re)

    # Step 8: x.re = x.re - t1
    x_re = fp_sub(x_re, t1)

    # Return Fp2 result
    return Fp2(re=x_re, im=x_im)

# Combine limbs (little endian → shift by 64 bits per limb)
def combine_limbs(limbs):
    val = 0
    for i, limb in enumerate(limbs):
        val += limb << (64 * i)
    return val

# Split big integer into 64-bit limbs (little endian)
def split_limbs(x, num_limbs=5):
    mask = (1 << 64) - 1
    limbs = []
    for i in range(num_limbs):
        limbs.append(x & mask)   # take lowest 64 bits
        x >>= 64                 # shift down for next limb
    return limbs



def combine_fp2(fp2):
    """Return combined integer values for an Fp2 element."""
    return {
        "re": combine_limbs(fp2.re),
        "im": combine_limbs(fp2.im),
    }

def combine_fp2_class(fp2):
    """Return combined integer values of an Fp2 class object."""
    res = combine_fp2(fp2)
    return Fp2(re=res["re"], im=res["im"])

def combine_theta_point(tp):
    """Return combined integer values for every coordinate in a ThetaPoint."""
    return {
        "x": combine_fp2(tp.x),
        "y": combine_fp2(tp.y),
        "z": combine_fp2(tp.z),
        "t": combine_fp2(tp.t),
    }

def combine_theta_point_class(tp):
    """Return combined integer values of an Fp2 class object."""
    res = {
        "x": combine_fp2_class(tp.x),
        "y": combine_fp2_class(tp.y),
        "z": combine_fp2_class(tp.z),
        "t": combine_fp2_class(tp.t),
    }
    return ThetaPoint(x=res["x"], y=res["y"], z=res["z"], t=res["t"])

def hadamard_sqisign(in1): # input is a ThetaPoint
    t1 = fp2_add(in1.x ,in1.y)
    t2 = fp2_sub(in1.x, in1.y)
    t3 = fp2_add(in1.z, in1.t)
    t4 = fp2_sub(in1.z, in1.t)

    out1 = ThetaPoint(x=0, y=0, z=0, t=0)
    out1.x = fp2_add(t1, t3)
    out1.y = fp2_add(t2, t4)
    out1.z = fp2_sub(t1, t3)
    out1.t = fp2_sub(t2, t4)

    return out1


def fp_sqr(x: Fp2)->Fp2:
    return fp_mul(x ,x)

def mont_decode(x_mont: int, q: int, R_inv: int) -> int:
    # x_mont is x*R mod q  -> return x
    return (x_mont * R_inv) % q

def mont_decode(x_mont: int, q: int, R_inv: int) -> int:
    # x_mont is x*R mod q  -> return x
    return (x_mont * R_inv) % q

def fp_encode_lsb(x: int, q: int, R_inv: int, is_mont: bool) -> int:
    """
    Return the exact bit that C reads as tmp_bytes[0] & 1 after fp_encode().
    """
    xn = mont_decode(x, q, R_inv) if is_mont else (x % q)
    return xn & 1




def fp_sqrt_canonical(z: int, q: int, R_inv: int, is_mont: bool) -> int:
    z %= q
    y = pow(z, (q + 1) // 4, q)   # candidate sqrt in normal arithmetic

    # If your field uses Montgomery internally, y here is normal already because pow() is normal.
    # But if you replace pow() with montgomery exponentiation later, keep this LSB logic.
    if fp_encode_lsb(y, q, R_inv, is_mont) == 1:
        y = (-y) % q
    return y


# def fp2_sqrt_match_c(a_in: Fp2, q: int, is_mont: bool = False) -> Fp2:
#     """
#     Matches the C fp2_sqrt() including:
#       - canonical fp_sqrt() output (even low bit)
#       - special-case when a1 == 0: force delta = a0
#       - final canonical sign rule on (re, im)
#     Assumes q is an odd prime with q % 4 == 3 and i^2 = -1.
#     """
#     a0 = a_in.re * R % q
#     a1 = a_in.im * R % q

#     a0 = a0 % q
#     a1 = a1 % q

#     # x0 = delta = sqrt(a0^2 + a1^2) in F_q, canonicalized like gf5248_sqrt
#     n = (a0 * a0 + a1 * a1) % q
#     # (Optional) If you want same assumptions as C fp_sqrt: it doesn't fail hard if non-square.
#     # Keep your residue check if you want safety.
#     # if n != 0 and pow(n, (q-1)//2, q) != 1:
#     #     raise ValueError("not a square")

#     delta = fp_sqrt_canonical(n, q, R_INV, is_mont)

#     # If a1 == 0, force delta = a0 (C does this to avoid delta = -a0 -> x0=0 later)
#     if a1 == 0:
#         delta = a0

#     # x0 = delta + a0, t0 = 2*x0
#     x0 = (delta + a0) % q
#     t0 = (2 * x0) % q

#     # x1 = t0^((q-3)/4)
#     x1 = pow(t0, (q - 3) // 4, q)

#     # x0 = x0*x1 ; x1 = x1*a1 ; t1 = (2*x0)^2
#     x0 = (x0 * x1) % q
#     x1 = (x1 * a1) % q
#     t1 = (2 * x0) % q
#     t1 = (t1 * t1) % q

#     # If t1 == t0 return (x0, x1) else return (x1, -x0)
#     if (t0 - t1) % q == 0:
#         re, im = x0, x1
#     else:
#         re, im = x1, (-x0) % q

#     # FINAL canonicalization exactly like C using fp_encode LSBs
#     re_is_zero = (re % q) == 0
#     re_is_odd  = fp_encode_lsb(re, q, R_INV, is_mont) == 1
#     im_is_odd  = fp_encode_lsb(im, q, R_INV, is_mont) == 1

#     negate_output = re_is_odd or (re_is_zero and im_is_odd)
#     if negate_output:
#         re = (-re) % q
#         im = (-im) % q

#     return Fp2(re, im)


def fp_select(a0, a1, ctl):
    return a0 if ctl == 0 else a1

def fp2_select(a0, a1, ctl):
    d = Fp2()
    d.re = fp_select(a0.re, a1.re, ctl)
    d.im = fp_select(a0.im, a1.im, ctl)

    return d




def fp2_exp3div4(a):
    """
    Implements the exact addition-chain shown:

      - computes intermediates like a^(2^k - 1)
      - ends with: a^(5*(2^246 - 1) + 4)

    Returns the final 'a' value (new FP element).
    """

    # fp_sqr(&z4, a); // a^2
    z4 = fp2_sqr(a)                      # a^2

    # fp_mul(&tmp, a, &z4); // a^3
    tmp = fp2_mul(a, z4)                 # a^3


    # fp_sqr(&z4, &z4); // a^4
    z4 = fp2_sqr(z4)                     # a^4

    # fp_mul(&t3, &tmp, &z4); // a^7
    t3 = fp2_mul(tmp, z4)                # a^7  = a^(2^3-1)


    # ---- Compute a^(2^6 - 1) = a^63 ----
    # fp_sqr(&t6, &t3); // a^14
    t6 = fp2_sqr(t3)                     # a^14
    # for (i=1; i<3; i++) fp_sqr(&t6,&t6);  (2 squarings)
    for _ in range(1, 3):
        t6 = fp2_sqr(t6)                 # a^28, a^56

    # fp2_mul(&t6, &t6, &t3); // a^63
    t6 = fp2_mul(t6, t3)                 # a^63 = a^(2^6-1)

    # ---- Compute a^(2^12 - 1) ----
    # fp2_sqr(a, &t6); then 5 more squarings (total 6)
    #t6 =  # copy t6 to a for squarings
    a =  Fp2(re=t6.re, im=t6.im) 
    for _ in range(0, 6):
        a = fp2_sqr(a)

    # fp2_mul(a, a, &t6)
    a = fp2_mul(a, t6)                   # a^(2^12-1)
    # ---- Compute a^(2^15 - 1) ----
    # 3 squarings
    for _ in range(0, 3):
        a = fp2_sqr(a)
    # * a^7
    a = fp2_mul(a, t3)                   # a^(2^15-1)
    # ---- Compute a^(2^30 - 1) ----
    tmp = fp2_sqr(a)
    for _ in range(1, 15):
        tmp = fp2_sqr(tmp)               # tmp = a^(2^30-2^15)
    a = fp2_mul(a, tmp)                  # a^(2^30-1)
    # ---- Compute a^(2^60 - 1) ----
    tmp = fp2_sqr(a)
    for _ in range(1, 30):
        tmp = fp2_sqr(tmp)
    a = fp2_mul(a, tmp)                  # a^(2^60-1)
    # ---- Compute a^(2^120 - 1) ----
    tmp = fp2_sqr(a)
    for _ in range(1, 60):
        tmp = fp2_sqr(tmp)
    a = fp2_mul(a, tmp)                  # a^(2^120-1)
    # ---- Compute a^(2^123 - 1) ----
    for _ in range(0, 3):
        a = fp2_sqr(a)
    a = fp2_mul(a, t3)                   # a^(2^123-1)
    # ---- Compute a^(2^246 - 1) ----
    tmp = fp2_sqr(a)
    for _ in range(1, 123):
        tmp = fp2_sqr(tmp)
    a = fp2_mul(a, tmp)                  # a^(2^246-1)
    
    # ---- Compute a^(5*(2^246 - 1)) ----
    tmp = fp2_sqr(a)                     # a^2
    tmp = fp2_sqr(tmp)                   # a^4
    a = fp2_mul(a, tmp)                  # a^5  (times previous exponent) => a^(5*(2^246-1))
    # ---- Compute a^(5*(2^246 - 1) + 4) ----
    a = fp2_mul(a, z4)                   # z4 is still a^4 from the start
    return a

# Assume FP_ENCODED_BYTES is defined in your codebase.
# FP2 is represented as a tuple (re, im)


def fp_sqrt_only_fp2(a_mont: Fp2):
    res = Fp2(ONE, 0)
    base = a_mont

    a_2 = fp2_mul(a_mont, a_mont)
    a_4 = fp2_mul(a_2, a_2)
    a_5 = fp2_mul(a_4, a_mont)

    res = a_5

    for i in range(246):
        res = fp2_mul(res, res)
    return res



def fp2_sqrt(a: Fp2):
    """
    Compute sqrt(a) in FP2.

    Parameters
    ----------
    a : tuple
        (re, im) in your FP representation.

    Returns
    -------
    tuple
        (re, im) = sqrt(a) in FP2, canonically signed like the C code.
    """
    a0 = Fp2(a.re, 0)
    a1 = Fp2(a.im, 0)  # re, im

    # x0 = δ = sqrt(a0^2 + a1^2)
    x0 = fp2_sqr(a0)
    x1 = fp2_sqr(a1)

    x0 = fp2_add(x0, x1)
    x0 = fp_sqrt_only_fp2(x0)

    # If a1 == 0, force δ = a0 (avoid δ = -a0 leading to x0=0 later)
    #print("fp_is_zero(a1.re, q): ", fp_is_zero(a1.re, q))
    x0 = Fp2(fp_select(x0.re, a0.re, fp_is_zero(a1.re, q)), 0)

    # x0 = δ + a0 ; t0 = 2*x0
    x0 = fp2_add(x0, a0)
    t0 = fp2_add(x0, x0)

    # x1 = t0^((p-3)/4)
    x1 = fp2_exp3div4(t0)

    # x0 = x0*x1 ; x1 = x1*a1 ; t1 = (2*x0)^2
    x0 = fp2_mul(x0, x1)
    x1 = fp2_mul(x1, a1)
    t1 = fp2_add(x0, x0)
    t1 = fp2_sqr(t1)

    # If t1 == t0, output is (x0, x1); else (x1, -x0)
    diff = fp2_sub(t0, t1)
    f = fp_is_zero(diff.re, q)

    neg_x0 = fp2_neg(x0)

    # t0_out = f ? x0 : x1
    # t1_out = f ? x1 : -x0
    t0_out = Fp2(fp_select(x1.re, x0.re, f), 0)
    t1_out = Fp2(fp_select(neg_x0.re, x1.re, f), 0)

    # Canonical sign fix:
    # negate if (t0 is odd) or (t0==0 and t1 is odd)
    t0_is_zero = fp_is_zero(t0_out.re, q)

    b0 = fp_encode(t0_out.re)  # bytes
    t0_is_odd = 1 if (b0[0] & 1) else 0

    b1 = fp_encode(t1_out.re)
    t1_is_odd = 1 if (b1[0] & 1) else 0

    negate_output = (t0_is_odd | (t0_is_zero & t1_is_odd))

    re = fp_select(t0_out.re, fp_neg(t0_out.re), negate_output)
    im = fp_select(t1_out.re, fp_neg(t1_out.re), negate_output)

    res = Fp2(re, im)

    return res




def fp2_sqr(y: Fp2) -> Fp2:
    """
    Computes x = y^2 in the quadratic extension field Fp2.

    Args:
        y (Fp2): input element (with fields re, im)
    Returns:
        Fp2: output squared element
    """

    # sum = y.re + y.im
    sum_ = fp_add(y.re, y.im)

    # diff = y.re - y.im
    diff = fp_sub(y.re, y.im)

    # x.im = (2* y.re) * y.im
    re_mul2 = fp_add(y.re, y.re)  # multiply by 2 via addition
    im_part = fp_mul(re_mul2, y.im)
   
    # x.re = (y.re + y.im) * (y.re - y.im)
    re_part = fp_mul(sum_, diff)

    # Return new Fp2 element
    return Fp2(re=re_part, im=im_part)



def fp2_sqr_c0(y: Fp2) -> Fp2:
    # sum = y.re + y.im
    sum_ = fp_add_sqr(y.re, y.im)

    # diff = y.re - y.im
    diff = fp_sub_sqr(y.re, y.im)

    # x.re = (y.re + y.im) * (y.re - y.im)
    re_part = fp_mul_fp_sqr(sum_, diff)

    # Return new Fp2 element
    return re_part

def fp2_sqr_c1(y: Fp2) -> Fp2:

    im_part = fp_add(y.re, y.re)  # multiply by 2 via addition
    im_part = fp_mul_fp_sqr(im_part, y.im)
    

    # Return new Fp2 element
    return im_part



def fp2_sqr_new(a: Fp2)->Fp2:
    re_part = fp2_sqr_c0(a)
    im_part = fp2_sqr_c1(a)

    return Fp2(re=re_part, im=im_part)


def pointwise_square(tp_in):
    out1 = ThetaPoint(x=0, y=0, z=0, t=0)
    out1.x = fp2_sqr(tp_in.x)
    out1.y = fp2_sqr(tp_in.y) 
    out1.z = fp2_sqr(tp_in.z)
    out1.t = fp2_sqr(tp_in.t)
    
    return out1


def to_squared_theta(input1):
    out11 = pointwise_square(input1)
    ret1 = hadamard_sqisign(out11)
    
    return ret1

def hadamard_and_ptwise_sq(input1):
    out11 = hadamard_sqisign(input1)
    ret1 = to_squared_theta(out11)

    return ret1

def theta_precomputation(A: ThetaStructure) -> ThetaStructure:
    """
    Python equivalent of:
        void theta_precomputation(theta_structure_t *A)

    Returns a *new* ThetaStructure A_new with computed fields (no in-place mutation).
    """

    # If already precomputed, just return as-is
    if A.precomputation:
        return A

    # --- Step 1: Compute A_dual = to_squared_theta(A.null_point)
    A_dual = to_squared_theta(A.null_point)

    # --- Step 2: Temporary fp2 intermediates
    t1 = fp2_mul(A_dual.x, A_dual.y)
    t2 = fp2_mul(A_dual.z, A_dual.t)

    # --- Step 3: Compute upper-case constants
    XYZ0 = fp2_mul(t1, A_dual.z)
    XYT0 = fp2_mul(t1, A_dual.t)
    YZT0 = fp2_mul(t2, A_dual.y)
    XZT0 = fp2_mul(t2, A_dual.x)

    # --- Step 4: Recompute t1, t2 using null_point
    t1 = fp2_mul(A.null_point.x, A.null_point.y)
    t2 = fp2_mul(A.null_point.z, A.null_point.t)

    # --- Step 5: Compute lower-case constants
    xyz0 = fp2_mul(t1, A.null_point.z)
    xyt0 = fp2_mul(t1, A.null_point.t)
    yzt0 = fp2_mul(t2, A.null_point.y)
    xzt0 = fp2_mul(t2, A.null_point.x)

    # --- Step 6: Build and return new structure
    A_new = ThetaStructure(
        null_point=A.null_point,
        precomputation=True
    )

    # Assign all computed constants
    A_new.XYZ0 = XYZ0
    A_new.XYT0 = XYT0
    A_new.YZT0 = YZT0
    A_new.XZT0 = XZT0
    A_new.xyz0 = xyz0
    A_new.xyt0 = xyt0
    A_new.yzt0 = yzt0
    A_new.xzt0 = xzt0

    return A_new








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



# Parameters / types
NWORDS_ORDER = 4                  # number of 64-bit limbs
RADIX = 64                        # limb size in bits
SQIsign_response_length = 126
HD_extra_torsion = 2
TORSION_EVEN_POWER = 248

BASIS_E0_PX = Fp2(re=0x024e4cc21a236db3f59a8d87ad37c0da8c8505452654b56d052b795624001810, im=0x026ef8a8622cda106a9bbbb814c83cbd5cc5efa9da1d4a82cd9d72c0cb907df8)
BASIS_E0_QX = Fp2(im=0x04ff74b2ac0249ecc63b049a9c3405e7e7b1cb210f2d30b26ad43baab72f065f, re=0x04f761f96b4a5f40e121cbd7b1571ed86634424982edefcc606e8b2029222fc7)

p_cofactor_for_2f = 5
P_COFACTOR_FOR_2F_BITLENGTH = 3

MASK64 = (1 << 64) - 1


@dataclass
class ECPoint:
    x: Fp2
    z: Fp2

@dataclass
class ECCurve:
    A: Fp2
    C: Fp2
    A24: Optional[ECPoint]
    is_A24_computed_and_normalized: bool

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
class PublicKey:
    curve: ECCurve
    hint_pk: int

@dataclass
class Signature:
    E_aux_A: Fp2
    mat_Bchall_can_to_B_chall: list   # 2x2 matrix of ints
    backtracking: int
    two_resp_length: int
    hint_aux: int
    hint_chall: int
    chall_coeff: int