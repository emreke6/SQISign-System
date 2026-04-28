from arith_ops import *
from basis import *
from ec import *

from our_model import *

import random

def generate_xmul_ladder(my_rom, start_idx, kbits=3):
    """
    Generates the Montgomery Ladder ROM instructions.
    Returns the next available index in the ROM.
    """
    i = start_idx
    
    # Initialize Loop
    my_rom[i] = pack(C_SET_LOOP, S_NONE, kbits); i += 1

    # --- LADDER LOOP START ---
    loop_start_pc = i
    my_rom[i] = pack(C_NEXT_BIT, S_NONE, 0); i += 1 

    # --- EXPLICIT CONDITIONAL SWAP (.x) ---
    my_rom[i] = pack(C_READ, S_REGA, 0); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 2); i += 1
    my_rom[i] = pack(C_CSWAP_REGS, S_NONE, 0); i += 1 
    my_rom[i] = pack(C_WD, S_REGA, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGB, 2); i += 1

    # --- EXPLICIT CONDITIONAL SWAP (.z) ---
    my_rom[i] = pack(C_READ, S_REGA, 1); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 3); i += 1
    my_rom[i] = pack(C_CSWAP_REGS, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGA, 1); i += 1
    my_rom[i] = pack(C_WD, S_REGB, 3); i += 1

    # --- xDBLADD CORE ---
    # t0 = P.x + P.z
    my_rom[i] = pack(C_READ, S_REGA, 0); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 1); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1

    # t1 = P.x - P.z
    my_rom[i] = pack(C_READ, S_REGA, 0); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 1); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 8); i += 1

    # Rx = t0^2
    my_rom[i] = pack(C_READ, S_REGA, 7); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 9); i += 1

    # t2 = Q.x - Q.z
    my_rom[i] = pack(C_READ, S_REGA, 2); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 3); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1


    # Sx = Q.x + Q.z
    my_rom[i] = pack(C_READ, S_REGA, 2); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 3); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 11); i += 1
    

    # t0 = t0 * t2
    my_rom[i] = pack(C_READ, S_REGA, 7); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 10); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1

    # Rz = t1^2
    my_rom[i] = pack(C_READ, S_REGA, 8); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 8); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 12); i += 1

    # t1 = t1 * Sx
    my_rom[i] = pack(C_READ, S_REGA, 8); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 11); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 8); i += 1

    # t2 = Rx - Rz
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 12); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1

    # R0.x = Rx * Rz
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 12); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 0); i += 1

    print("Generated xDBLADD core instructions up to Rx computation 44. ", i)

    # Sx2 = A24.x * t2
    my_rom[i] = pack(C_READ, S_REGA, 4); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 10); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 13); i += 1

    print("Generated xDBLADD core instructions up to Rx computation 33. ", i)

    # Sz = t0 - t1
    my_rom[i] = pack(C_READ, S_REGA, 7); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 8); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 14); i += 1

    # Rz = Rz + Sx2
    my_rom[i] = pack(C_READ, S_REGA, 12); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 13); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 12); i += 1

    # Sx = t0 + t1
    my_rom[i] = pack(C_READ, S_REGA, 7); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 8); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 11); i += 1

    # R0.z = Rz * t2
    my_rom[i] = pack(C_READ, S_REGA, 12); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 10); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 1); i += 1

    # Sz = Sz^2
    my_rom[i] = pack(C_READ, S_REGA, 14); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 14); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 14); i += 1

    # Sx = Sx^2
    my_rom[i] = pack(C_READ, S_REGA, 11); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 11); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 11); i += 1

    # R1.z = Sz * PQ.x
    my_rom[i] = pack(C_READ, S_REGA, 14); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 5); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 3); i += 1

    # R1.x = Sx * PQ.z
    my_rom[i] = pack(C_READ, S_REGA, 11); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGC, 2); i += 1

    # Loop Back
    my_rom[i] = pack(C_DEC_JNZ, S_NONE, loop_start_pc); i += 1

    # --- FINAL EXPLICIT SWAP CORRECTION ---
    my_rom[i] = pack(C_FINAL_SWAP, S_NONE, 0); i += 1

    # Final Swap (.x)
    my_rom[i] = pack(C_READ, S_REGA, 0); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 2); i += 1
    my_rom[i] = pack(C_CSWAP_REGS, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGA, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGB, 2); i += 1

    # Final Swap (.z)
    my_rom[i] = pack(C_READ, S_REGA, 1); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 3); i += 1
    my_rom[i] = pack(C_CSWAP_REGS, S_NONE, 0); i += 1
    my_rom[i] = pack(C_WD, S_REGA, 1); i += 1
    my_rom[i] = pack(C_WD, S_REGB, 3); i += 1

    return i # Return the next PC address


def generate_difference_point_op(my_rom, PC_in):
    # MEM-MAP:
    # P.x = 0
    # P.z = 1
    # Q.x = 2
    # Q.z = 3
    # curve_A = 4
    # curve_C = 5 
    # t0 = 6
    # t1 = 7
    # Bxx = 8
    # Bxz = 9
    # Bzz = 10

    # Temporary variables
    # t0 = fp2_mul(P.x, Q.x)
    i = PC_in
    my_rom[i] = pack(C_READ, S_REGA, 0); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 2); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1
    
    #t1 = fp2_mul(P.z, Q.z)
    my_rom[i] = pack(C_READ, S_REGA, 1); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 3); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1


    # Bxx = fp2_sub(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    # Bxx is now in 8
    my_rom[i] = pack(C_WD, S_REGC, 8); i += 1

    # Bxz = fp2_add(t0, t1)
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    # Bxz is now in 9
    my_rom[i] = pack(C_WD, S_REGC, 9); i += 1

    # Bxx = fp2_sqr(Bxx)
    my_rom[i] = pack(C_READ, S_REGA, 8); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 8); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bxx is now in 8
    my_rom[i] = pack(C_WD, S_REGC, 8); i += 1


    #Bxx = fp2_mul(Bxx, curve.C)
    my_rom[i] = pack(C_READ, S_REGA, 8); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 5); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bxx is now in 8
    my_rom[i] = pack(C_WD, S_REGC, 8); i += 1

    

    #t0 = fp2_mul(P.x, Q.z)
    my_rom[i] = pack(C_READ, S_REGA, 0); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 3); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1
    
    # t1 = fp2_mul(P.z, Q.x)
    my_rom[i] = pack(C_READ, S_REGA, 1); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 2); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1

    #Bzz = fp2_add(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    # Bzz is now in 10
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1

    
    #Bxz = fp2_mul(Bxz, Bzz)
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 10); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bxz is now in 9
    my_rom[i] = pack(C_WD, S_REGC, 9); i += 1


    #Bzz = fp2_sub(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    # Bzz is now in 10
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1

    
    #Bzz = fp2_sqr(Bzz)
    my_rom[i] = pack(C_READ, S_REGA, 10); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 10); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bzz is now in 10
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1
    
    
    #Bzz = fp2_mul(Bzz, curve.C)
    my_rom[i] = pack(C_READ, S_REGA, 10); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 5); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bzz is now in 10
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1

    #Bxz = fp2_mul(Bxz, curve.C)
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 5); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bxz is now in 9
    my_rom[i] = pack(C_WD, S_REGC, 9); i += 1

    #t0 = fp2_mul(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1
    
    #t0 = fp2_mul(t0, curve.A)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 4); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    #t0 = fp2_add(t0, t0)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1
    

    #Bxz = fp2_add(Bxz, t0)
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_ADD, S_NONE, 0); i += 1
    # Bxz is now in 9
    my_rom[i] = pack(C_WD, S_REGC, 9); i += 1

    

    # Normalization factor
    #t0 = Fp2(curve.C.re, fp_neg(curve.C.im))
    my_rom[i] = pack(C_READ, S_REGA, 5); i += 1
    my_rom[i] = pack(C_CONJ, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    # t0 = fp2_sqr(t0)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1


    #t0 = fp2_mul(t0, curve.C)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 5); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    #t1 = Fp2(P.z.re, fp_neg(P.z.im))
    my_rom[i] = pack(C_READ, S_REGA, 1); i += 1
    my_rom[i] = pack(C_CONJ, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1


    #t1 = fp2_sqr(t1)
    my_rom[i] = pack(C_READ, S_REGA, 7); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1


    #t0 = fp2_mul(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    #t1 = Fp2(Q.z.re, fp_neg(Q.z.im))
    my_rom[i] = pack(C_READ, S_REGA, 3); i += 1
    my_rom[i] = pack(C_CONJ, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1


    #t1 = fp2_sqr(t1)
    my_rom[i] = pack(C_READ, S_REGA, 7); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1
    
    
    #t0 = fp2_mul(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    #Bxx = fp2_mul(Bxx, t0)
    my_rom[i] = pack(C_READ, S_REGA, 8); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bxx is now in 8
    my_rom[i] = pack(C_WD, S_REGC, 8); i += 1
    
    #Bxz = fp2_mul(Bxz, t0)
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bxz is now in 9
    my_rom[i] = pack(C_WD, S_REGC, 9); i += 1
    
    
    #Bzz = fp2_mul(Bzz, t0)
    my_rom[i] = pack(C_READ, S_REGA, 10); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 6); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # Bzz is now in 10
    my_rom[i] = pack(C_WD, S_REGC, 10); i += 1

    # Solve quadratic
    #t0 = fp2_sqr(Bxz)
    my_rom[i] = pack(C_READ, S_REGA, 9); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 9); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    
    
    #t1 = fp2_mul(Bxx, Bzz)
    my_rom[i] = pack(C_READ, S_REGA, 8); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 10); i += 1
    my_rom[i] = pack(C_MUL, S_NONE, 0); i += 1
    # t1 is now in 7
    my_rom[i] = pack(C_WD, S_REGC, 7); i += 1


    #t0 = fp2_sub(t0, t1)
    my_rom[i] = pack(C_READ, S_REGA, 6); i += 1
    my_rom[i] = pack(C_READ, S_REGB, 7); i += 1
    my_rom[i] = pack(C_SUB, S_NONE, 0); i += 1
    # t0 is now in 6
    my_rom[i] = pack(C_WD, S_REGC, 6); i += 1

    return i # Return the next PC address after the generated code


def random_ec_point():
    random_point = ec_point_init()

    random_point.x = Fp2(random.randint(0, q-1), random.randint(0, q-1))
    random_point.z = Fp2(random.randint(0, q-1), random.randint(0, q-1))

    return random_point

def random_ec_basis():
    random_basis = ECBasis(P=ec_point_init(), Q=ec_point_init(), PmQ=ec_point_init())

    return random_basis


def random_ec_curve():
    random_curve = ECCurve(A=Fp2(0, 0), C=Fp2(0, 0), A24=Fp2(0, 0), is_A24_computed_and_normalized=False)

    random_curve.A = Fp2(random.randint(0, q-1), random.randint(0, q-1))
    random_curve.C = Fp2(random.randint(0, q-1), random.randint(0, q-1))
    random_curve.A24 = random_ec_point()
    random_curve.is_A24_computed_and_normalized = True

    return random_curve

if __name__ == "__main__":

    random.seed(42)  # For reproducibility

    PQ2 = random_ec_basis()
    PQ2.P.z = fp2_set_one()  # Ensure P.z is set to 1 for valid projective coordinates
    curve_test = random_ec_curve()
    curve_test.is_A24_computed_and_normalized = True
    f = random.randint(0, q-1)
    hint = random.randint(0, 100)

    # Test the basis conversion function
    #ret, curve_new, PQ2_res = ec_curve_to_basis_2f_from_hint(PQ2, curve_test, f, hint)


    # print("Test result:", ret)
    # print_fp2("New curve A:", curve_new.A)
    # print_fp2("New curve C:", curve_new.C)
    # print(curve_new)
    # print_hex_point("New curve A24:", curve_new.A24)
    # print_fp2("PQ2.P.x:", PQ2_res.P.x)
    # print_fp2("PQ2.P.z:", PQ2_res.P.z)
    # print_fp2("PQ2.Q.x:", PQ2_res.Q.x)
    # print_fp2("PQ2.Q.z:", PQ2_res.Q.z)
    # print_fp2("PQ2.PmQ.x:", PQ2_res.PmQ.x)
    # print_fp2("PQ2.PmQ.z:", PQ2_res.PmQ.z)


    memory_address_map_ec_curve_basis = {
        "curve_A": 64,
        "curve_C": 65,
        "curve_A24_x": 66,
        "PQ2_P_x": 67,
        "PQ2_P_z": 68,
        "PQ2_Q_x": 69,
        "PQ2_Q_z": 70,
        "PQ2_PmQ_x": 71,
        "PQ2_PmQ_z": 72,
        "hint_P": 75,
        "ONE": 128,
        "ZERO": 129,
    }

    my_rom = [pack(C_HALT, S_NONE, 0)] * 4096  # Initialize ROM with HALT instructions


    # ================================================================
    # STEP 1: OFFLOAD FROM GLOBAL TO LOCAL
    # ================================================================
    # We bring CURVE_A into a local scratchpad (Addr 0)

    my_rom[0]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["curve_A"]); my_rom[1]  = pack(C_WD, S_REGA, 0) # load curve_a
    my_rom[2]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["hint_P"]); my_rom[3]  = pack(C_WD, S_REGA, 1) # load hint_P
    my_rom[4]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_x"]); my_rom[5]  = pack(C_WD, S_REGA, 2) # load Px
    my_rom[6]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_z"]); my_rom[7]  = pack(C_WD, S_REGA, 3) # load Pz

    # fp_mul(curve_a, hint_P) -> regC
    my_rom[8]  = pack(C_READ, S_REGA, 0)
    my_rom[9]  = pack(C_READ, S_REGB, 1)
    my_rom[10] = pack(C_MUL, S_NONE, 0) 
    my_rom[11] = pack(C_WD, S_REGC, memory_address_map_ec_curve_basis["PQ2_P_x"]) # fp_mul(curve_a, hint_P) -> regC

    # Q.x = fp2_add(curve.A, P.x)
    my_rom[12] = pack(C_READ, S_REGA, 0)
    my_rom[13] = pack(C_READ, S_REGB, 2)
    my_rom[14] = pack(C_ADD, S_NONE, 0)
    my_rom[15] = pack(C_WD, S_REGC, memory_address_map_ec_curve_basis["PQ2_Q_x"]) # fp2_add(curve.A, P.x) -> Q.x

    # Q.x = fp2_neg(Q.x)
    my_rom[16] = pack(C_RESET_REGS, S_REG_A_B, 0)
    my_rom[17] = pack(C_READ, S_REGB, memory_address_map_ec_curve_basis["PQ2_Q_x"])
    my_rom[18] = pack(C_SUB, S_NONE, 0)
    my_rom[19] = pack(C_WD, S_REGC, memory_address_map_ec_curve_basis["PQ2_Q_x"])

    my_rom[20] = pack(C_SET_ONE_REGA, S_REG_A_B, 0) # fp2_set_one()
    my_rom[21] = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_z"]) # fp2_set_one()


    # PREAMBLE FOR XMUL OF P
    my_rom[22]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["ONE"]); my_rom[23]  = pack(C_WD, S_REGA, 0) # load R0
    my_rom[24]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["ZERO"]); my_rom[25]  = pack(C_WD, S_REGA, 1) # load R1
    my_rom[26]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_x"]); my_rom[27]  = pack(C_WD, S_REGA, 2) # load P.x
    my_rom[28]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_z"]); my_rom[29]  = pack(C_WD, S_REGA, 3) # load P.z
    my_rom[30]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["curve_A24_x"]); my_rom[31]  = pack(C_WD, S_REGA, 4) # load A24.x
    my_rom[32]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_x"]); my_rom[33]  = pack(C_WD, S_REGA, 5) # load P.x
    my_rom[34]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_z"]); my_rom[35]  = pack(C_WD, S_REGA, 6) # load P.z

    PC = generate_xmul_ladder(my_rom, 36, kbits=3)

    my_rom[PC]  = pack(C_READ, S_REGA, 0); my_rom[PC+1]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_x"]) # load R0
    my_rom[PC+2]  = pack(C_READ, S_REGA, 1); my_rom[PC+3]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_z"]) # load R1
    PC = PC + 4


    # PREAMBLE FOR XMUL OF Q
    my_rom[PC]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["ONE"]); my_rom[PC+1]  = pack(C_WD, S_REGA, 0) # load R0
    my_rom[PC+2]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["ZERO"]); my_rom[PC+3]  = pack(C_WD, S_REGA, 1) # load R1
    my_rom[PC+4]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_x"]); my_rom[PC+5]  = pack(C_WD, S_REGA, 2) # load Q.x
    my_rom[PC+6]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_z"]); my_rom[PC+7]  = pack(C_WD, S_REGA, 3) # load Q.z
    my_rom[PC+8]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["curve_A24_x"]); my_rom[PC+9]  = pack(C_WD, S_REGA, 4) # load A24.x
    my_rom[PC+10]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_x"]); my_rom[PC+11]  = pack(C_WD, S_REGA, 5) # load Q.x
    my_rom[PC+12]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_z"]); my_rom[PC+13]  = pack(C_WD, S_REGA, 6) # load Q.z

    PC = PC + 14

    PC = generate_xmul_ladder(my_rom, PC, kbits=3)

    # POSTAMBLE
    my_rom[PC]  = pack(C_READ, S_REGA, 0); my_rom[PC+1]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_x"]) # load R0
    my_rom[PC+2]  = pack(C_READ, S_REGA, 1); my_rom[PC+3]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_z"]) # load R1

    PC = PC + 4


    # PREAMBLE FOR DIFFERENCE POINT RECOVERY
    my_rom[PC]    = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_x"])
    my_rom[PC+1]  = pack(C_WD, S_REGA, 0)  # local P.x

    my_rom[PC+2]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_P_z"])
    my_rom[PC+3]  = pack(C_WD, S_REGA, 1)  # local P.z

    my_rom[PC+4]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_x"])
    my_rom[PC+5]  = pack(C_WD, S_REGA, 2)  # local Q.x

    my_rom[PC+6]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_z"])
    my_rom[PC+7]  = pack(C_WD, S_REGA, 3)  # local Q.z

    my_rom[PC+8]  = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["curve_A"])
    my_rom[PC+9]  = pack(C_WD, S_REGA, 4)  # local curve.A

    my_rom[PC+10] = pack(C_READ, S_REGA, memory_address_map_ec_curve_basis["curve_C"])
    my_rom[PC+11] = pack(C_WD, S_REGA, 5)  # local curve.C

    
    PC = PC + 12

    PC = generate_difference_point_op(my_rom, PC)

    # POSTAMBLE FOR DIFFERENCE POINT RECOVERY (CHANGE AFTERWARDS)
    my_rom[PC]  = pack(C_READ, S_REGA, 6); my_rom[PC+1]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_x"]) # load R0
    my_rom[PC+2]  = pack(C_READ, S_REGA, 7); my_rom[PC+3]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_Q_z"]) # load R1

    my_rom[PC+4]  = pack(C_READ, S_REGA, 6); my_rom[PC+5]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_PmQ_x"]) # load R0
    my_rom[PC+6]  = pack(C_READ, S_REGA, 7); my_rom[PC+7]  = pack(C_WD, S_REGA, memory_address_map_ec_curve_basis["PQ2_PmQ_z"]) # load R0

    PC = PC + 12




    
    my_rom[PC] = pack(C_HALT, S_NONE, 0) ; PC = PC+1
    




    # 1. Initialize Memory with Test Constants
    initial_mem = {
        memory_address_map_ec_curve_basis["curve_A"]: curve_test.A,    # Curve.A
        memory_address_map_ec_curve_basis["curve_C"]: curve_test.C,    # Curve.C
        memory_address_map_ec_curve_basis["curve_A24_x"]: curve_test.A24.x,    # Curve.A24.x
        memory_address_map_ec_curve_basis["PQ2_P_x"]: PQ2.P.x,    # 
        memory_address_map_ec_curve_basis["PQ2_P_z"]: PQ2.P.z,    # 
        memory_address_map_ec_curve_basis["PQ2_Q_x"]: PQ2.Q.x,    # 
        memory_address_map_ec_curve_basis["PQ2_Q_z"]: PQ2.Q.z,    # 
        memory_address_map_ec_curve_basis["PQ2_PmQ_x"]: PQ2.PmQ.x,    # 
        memory_address_map_ec_curve_basis["PQ2_PmQ_z"]: PQ2.PmQ.z,    # 
        memory_address_map_ec_curve_basis["hint_P"]: Fp2((hint>>1) * R % q, 0),    # hint_P
        memory_address_map_ec_curve_basis["ONE"]: Fp2(ONE, 0),    # ONE
        memory_address_map_ec_curve_basis["ZERO"]: Fp2(0, 0),    # ZERO
    }

    # 2. Run Engine
    cu = ControlUnitSim(initial_mem)
    engine = RomEngineSim(cu, my_rom, scalar=p_cofactor_for_2f)
    engine.run(verbose=True)

    # 3. Test Results
    print("\n--- RESULTS VALIDATION ---")
    px = cu.memory.get(memory_address_map_ec_curve_basis["PQ2_P_x"])

    res = fp2_mul_small(curve_test.A, hint>>1)
    res2 = fp2_add(curve_test.A, PQ2.P.x)
    res2 = fp2_neg(res2)




    print_fp2("Expected P.x:", res)
    print_fp2("Actual P.x:", px)

    

    print_fp2("Expected aa:", res2)
    print_fp2("Actual aa:", cu.memory.get(memory_address_map_ec_curve_basis["PQ2_Q_x"]))

    p_in = ECPoint(x=res, z=fp2_set_one())

    P_res = xMUL(p_in, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, curve_test)

    print_fp2("Expected aa2:", P_res.x)
    print_fp2("Actual aa2:", cu.memory.get(memory_address_map_ec_curve_basis["PQ2_P_x"]))

    q_in = ECPoint(x=res2, z=fp2_set_one())

    Q_res = xMUL(q_in, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, curve_test)

    print_fp2("Expected aa2:", Q_res.x)
    print_fp2("Actual aa2:", cu.memory.get(memory_address_map_ec_curve_basis["PQ2_Q_x"]))




    difference_point_res = difference_point(P_res, Q_res, curve_test)

    print_fp2("Expected difference point x:", difference_point_res.x)
    print_fp2("Actual difference point x:", cu.memory.get(6))

    print_fp2("Expected difference point x:", difference_point_res.z)
    print_fp2("Actual difference point x:", cu.memory.get(7))
