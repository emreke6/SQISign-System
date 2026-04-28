from arith_ops import *
from ec import *
from basis import *
from Crypto.Hash import SHAKE256



def ec_j_inv(curve: ECCurve):
    # j-invariant computation for Montgommery coefficient A2=(A+2C:4C)

    t1 = fp2_sqr(curve.C)
    j_inv = fp2_sqr(curve.A)
    t0 = fp2_add(t1, t1)
    t0 = fp2_sub(j_inv, t0)
    t0 = fp2_sub(t0, t1)
    j_inv = fp2_sub(t0, t1)
    t1 = fp2_sqr(t1)
    j_inv = fp2_mul(j_inv, t1)
    t0 = fp2_add(t0, t0)
    t0 = fp2_add(t0, t0)
    t1 = fp2_sqr(t0)
    t0 = fp2_mul(t0, t1)
    t0 = fp2_add(t0, t0)
    t0 = fp2_add(t0, t0)
    j_inv = fp2_inv(j_inv)
    j_inv = fp2_mul(t0, j_inv)

    return j_inv


import sha3


# -----------------------------------------------------------
# Incremental SHAKE256 wrapper (TweetFIPS202 compatible)
# -----------------------------------------------------------

class Shake256Inc:
    """
    Incremental SHAKE256 wrapper compatible with:
        keccak_inc_init / absorb / finalize / squeeze
    from TweetFIPS202-style C implementations.

    pysha3 already uses correct SHAKE domain separation (0x1F),
    so finalize() is just a logical barrier.
    """

    def __init__(self):
        self._ctx = sha3.shake_256()
        self._finalized = False

    def absorb(self, data: bytes):
        if self._finalized:
            raise RuntimeError("Cannot absorb after finalize()")
        if not isinstance(data, (bytes, bytearray, memoryview)):
            raise TypeError("absorb() needs bytes-like object")

        self._ctx.update(bytes(data))

    def finalize(self):
        # TweetFIPS finalize inserts padding;
        # pysha3 already handles it internally.
        self._finalized = True

    def squeeze(self, outlen: int) -> bytes:
        if not self._finalized:
            self.finalize()

        return self._ctx.digest(outlen)


# -----------------------------------------------------------
# Masking helper (matching mp_mod_2exp behavior)
# -----------------------------------------------------------

def mask_scalar_top_bits(buf: bytearray, nbits: int, radix=None):
    """
    Keep only lowest nbits bits (little-endian).
    Equivalent to mp_mod_2exp in most PQC reference code.
    """

    if nbits <= 0:
        return bytearray(len(buf))

    x = int.from_bytes(buf, "little")
    x &= (1 << nbits) - 1
    return bytearray(x.to_bytes(len(buf), "little"))


# -----------------------------------------------------------
# Main hash_to_challenge (your original logic)
# -----------------------------------------------------------

def hash_to_challenge(
    pk,
    com_curve,
    message: bytes,
    length: int,
    FP2_ENCODED_BYTES,
    SECURITY_BITS,
    HASH_ITERATIONS,
    TORSION_EVEN_POWER,
    SQIsign_response_length,
    RADIX,
):

    # --------------------------------------------------------
    # Step 1: j-invariant encoding
    # --------------------------------------------------------
    j1 = ec_j_inv(pk.curve)
    j2 = ec_j_inv(com_curve)

    b1 = fp2_encode(j1)
    b2 = fp2_encode(j2)

    assert len(b1) == FP2_ENCODED_BYTES
    assert len(b2) == FP2_ENCODED_BYTES

    msg = bytes(message[:length])

    # --------------------------------------------------------
    # First SHAKE
    # --------------------------------------------------------
    hash_bits = 2 * SECURITY_BITS
    hash_bytes = (hash_bits + 7) // 8

    ctx = Shake256Inc()
    ctx.absorb(b1)
    ctx.absorb(b2)
    ctx.absorb(msg)
    ctx.finalize()

    scalar = bytearray(ctx.squeeze(hash_bytes))
    scalar = mask_scalar_top_bits(scalar, hash_bits, RADIX)

    # --------------------------------------------------------
    # Iterated hashing
    # --------------------------------------------------------
    for _ in range(2, HASH_ITERATIONS):
        ctx = Shake256Inc()
        ctx.absorb(bytes(scalar))
        ctx.finalize()

        scalar = bytearray(ctx.squeeze(hash_bytes))
        scalar = mask_scalar_top_bits(scalar, hash_bits, RADIX)

    # --------------------------------------------------------
    # Final squeeze
    # --------------------------------------------------------
    final_bits = TORSION_EVEN_POWER - SQIsign_response_length
    final_bytes = (final_bits + 7) // 8

    ctx = Shake256Inc()
    ctx.absorb(bytes(scalar))
    ctx.finalize()

    out = bytearray(ctx.squeeze(final_bytes))
    out = mask_scalar_top_bits(out, final_bits, RADIX)

    # mp_mod_2exp equivalent
    out_int = int.from_bytes(out, "little") & ((1 << SECURITY_BITS) - 1)

    return out_int

