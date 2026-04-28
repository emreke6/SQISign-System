from arith_ops import *

class ControlUnitSim:
    def __init__(self, memory_dict=None):
        self.a = Fp2(re=0, im=0)  # a1, b1
        self.b = Fp2(re=0, im=0)  # a2, b2
        self.c = Fp2(re=0, im=0)  # c1, c2
        self.memory = memory_dict if memory_dict is not None else {}
        self.busy = False

    def step(self, op_code, data_select, data1, swap_flag=False):
        self.busy = True
        
        # Decode Command [8:4]
        cmd = (op_code >> 4) & 0x1F
        
        # CSWAP Logic (Triggered by C_CSWAP_REGS or op_code 0x08)
        if (op_code & 0x008) or (cmd == 0): # Handles the C_CSWAP_REGS case
            if swap_flag:
                self.a, self.b = self.b, self.a

        # Memory/Arithmetic Operations
        if cmd == 0b00001:   # RD
            val = self.memory.get(data1, Fp2(re=0, im=0))
            if data_select == 1: self.a = val
            elif data_select == 2: self.b = val
        elif cmd == 0b00111: # CONJ (New Command for Conjugation)
            # Conjugates Reg A and stores in Reg C
            self.c = Fp2(re=self.a.re, im=fp_neg(self.a.im))
        elif cmd == 0b00101: # RESET (New Command)
            self.a = Fp2(re=0, im=0)
            self.b = Fp2(re=0, im=0)
        elif cmd == 0b00110: # SET_ONE_REGA (New Command)
            self.a = Fp2(re=ONE, im=0)
        elif cmd == 0b00010: # ADD
            self.c = fp2_add(self.a, self.b)
        elif cmd == 0b00100: # SUB
            self.c = fp2_sub(self.a, self.b)
        elif cmd == 0b01000: # MUL
            self.c = fp2_mul(self.a, self.b)
        elif cmd == 0b10000: # WD
            if data_select == 1: self.memory[data1] = self.a
            elif data_select == 2: self.memory[data1] = self.b
            elif data_select == 3: self.memory[data1] = self.c

        self.busy = False

class RomEngineSim:
    def __init__(self, cu, rom_data, scalar=0b101):
        self.cu = cu
        self.rom = rom_data
        self.pc = 0
        self.loop_cnt = 0
        #self.scalar = 0b101  # Matches Verilog test scalar
        self.scalar_init = scalar
        self.scalar = scalar
        self.prev_bit = 0
        self.swap_reg = 0
        self.done = False

        # Opcode Constants (Matches Verilog localparams)
        self.C_HALT       = 0b00000_0000
        self.C_SET_LOOP   = 0b00000_0001
        self.C_DEC_JNZ    = 0b00000_0010
        self.C_NEXT_BIT   = 0b00000_0011
        self.C_FINAL_SWAP = 0b00000_0100
        self.C_CSWAP_REGS = 0b00000_1000
        self.C_RESET_REGS = 0b00101_0000
        self.C_SET_ONE_REGA = 0b00110_0000
        self.C_CONJ         = 0b00111_0000

    def reset_ladder_state(self, kbits):
        self.loop_cnt = kbits
        self.scalar = self.scalar_init
        self.prev_bit = 0
        self.swap_reg = 0
        
    def run(self, verbose=True):
        if verbose:
            print(f"{'PC':<4} | {'Instruction Type':<15} | {'Swap':<5} | {'Loop':<4} | {'Registers'}")
            print("-" * 70)

        while not self.done:
            instr = self.rom[self.pc]
            print(f"Executing instruction at PC={self.pc}: {bin(instr)}")  # Debug: Show raw instruction
            print("Decoded Instruction - ", end="")
            op_code = (instr >> 10) & 0x1FF
            data_sel = (instr >> 8) & 0x3
            data1 = instr & 0xFF

            print(f"Opcode: {bin(op_code)}, Data Select: {data_sel}, Data1: {data1}")  # Debug: Show decoded fields

            # --- Control Flow Logic ---
            if op_code == self.C_HALT:
                self.done = True
                type_name = "HALT"

            elif op_code == self.C_SET_LOOP:
                self.reset_ladder_state(data1)
                self.loop_cnt = data1
                self.pc += 1
                type_name = "SET_LOOP"

            elif op_code == self.C_DEC_JNZ:
                if self.loop_cnt > 1:
                    self.loop_cnt -= 1
                    self.pc = data1
                else:
                    self.pc += 1
                type_name = "DEC_JNZ"

            elif op_code == self.C_NEXT_BIT:
                # current_bit = (self.scalar >> 2) & 1
                # self.swap_reg = current_bit ^ self.prev_bit
                # self.prev_bit = current_bit
                # self.scalar = (self.scalar << 1) & 0x7
                current_bit = (self.scalar >> (self.loop_cnt - 1)) & 1

                self.swap_reg = current_bit ^ self.prev_bit
                self.prev_bit = current_bit
                self.pc += 1
                type_name = "NEXT_BIT"

            elif op_code == self.C_FINAL_SWAP:
                self.swap_reg = self.prev_bit
                self.pc += 1
                type_name = "FINAL_SWAP"

            # --- Datapath Logic ---
            else:
                self.cu.step(op_code, data_sel, data1, swap_flag=bool(self.swap_reg))
                self.pc += 1
                type_name = "DATAPATH"

            if verbose:
                reg_info = f"A:{self.cu.a.re}, B:{self.cu.b.re}, C:{self.cu.c.re}"
                #print(f"{self.pc-1:<4} | {type_name:<15} | {self.swap_reg:<5} | {self.loop_cnt:<4} | {reg_info}")

# --- Example Usage ---

# 1. Build the ROM using the packing logic from Verilog {Op(9), Sel(2), Addr(8)}
def pack(op, sel, addr):
    return (op << 10) | (sel << 8) | addr

# Recreating a snippet of your xDBLADD ROM
C_READ = 0b00001_0000; C_ADD = 0b00010_0000; C_WD = 0b10000_0000; C_HALT = 0

C_RESET_REGS = 0b00101_0000
C_SET_ONE_REGA = 0b00110_0000

# Constants for Packing
C_READ = 0b00001_0000
C_ADD  = 0b00010_0000
C_SUB  = 0b00100_0000
C_MUL  = 0b01000_0000
C_WD   = 0b10000_0000
C_HALT = 0b00000_0000

C_SET_LOOP   = 0b00000_0001
C_DEC_JNZ    = 0b00000_0010
C_NEXT_BIT   = 0b00000_0011
C_FINAL_SWAP = 0b00000_0100
C_CSWAP_REGS = 0b00000_1000
C_CONJ         = 0b00111_0000

S_NONE = 0
S_REGA = 1
S_REGB = 2
S_REGC = 3

S_REG_A_B = 4

def pack(op, sel, addr):
    return (op << 10) | (sel << 8) | addr

# Initialize 128-word ROM
my_rom = [pack(C_HALT, S_NONE, 0)] * 128

# --- Initial Data Loading (PC 0-13) ---
# Moving data from external memory (64+) to working scratchpad (0+)
my_rom[0]  = pack(C_READ, S_REGA, 64); my_rom[1]  = pack(C_WD, S_REGA, 0)
my_rom[2]  = pack(C_READ, S_REGA, 65); my_rom[3]  = pack(C_WD, S_REGA, 1)
my_rom[4]  = pack(C_READ, S_REGA, 66); my_rom[5]  = pack(C_WD, S_REGA, 2)
my_rom[6]  = pack(C_READ, S_REGA, 67); my_rom[7]  = pack(C_WD, S_REGA, 3)
my_rom[8]  = pack(C_READ, S_REGA, 68); my_rom[9]  = pack(C_WD, S_REGA, 4)
my_rom[10] = pack(C_READ, S_REGA, 69); my_rom[11] = pack(C_WD, S_REGA, 5)
my_rom[12] = pack(C_READ, S_REGA, 70); my_rom[13] = pack(C_WD, S_REGA, 6)

# --- xDBLADD ALGORITHM CORE ---

# Step 1: t0 = P.x + P.z
my_rom[14] = pack(C_READ, S_REGA, 0); my_rom[15] = pack(C_READ, S_REGB, 1)
my_rom[16] = pack(C_ADD,  S_NONE, 0); my_rom[17] = pack(C_WD,   S_REGC, 7)

# Step 2: t1 = P.x - P.z
my_rom[18] = pack(C_READ, S_REGA, 0); my_rom[19] = pack(C_READ, S_REGB, 1)
my_rom[20] = pack(C_SUB,  S_NONE, 0); my_rom[21] = pack(C_WD,   S_REGC, 8)

# Step 3: Rx = t0^2
my_rom[22] = pack(C_READ, S_REGA, 7); my_rom[23] = pack(C_READ, S_REGB, 7)
my_rom[24] = pack(C_MUL,  S_NONE, 0); my_rom[25] = pack(C_WD,   S_REGC, 9)

# Step 4: t2 = Q.x - Q.z
my_rom[26] = pack(C_READ, S_REGA, 2); my_rom[27] = pack(C_READ, S_REGB, 3)
my_rom[28] = pack(C_SUB,  S_NONE, 0); my_rom[29] = pack(C_WD,   S_REGC, 10)

# Step 5: Sx = Q.x + Q.z
my_rom[30] = pack(C_READ, S_REGA, 2); my_rom[31] = pack(C_READ, S_REGB, 3)
my_rom[32] = pack(C_ADD,  S_NONE, 0); my_rom[33] = pack(C_WD,   S_REGC, 11)

# Step 6: t0 = t0 * t2
my_rom[34] = pack(C_READ, S_REGA, 7); my_rom[35] = pack(C_READ, S_REGB, 10)
my_rom[36] = pack(C_MUL,  S_NONE, 0); my_rom[37] = pack(C_WD,   S_REGC, 7)

# Step 7: Rz = t1^2
my_rom[38] = pack(C_READ, S_REGA, 8); my_rom[39] = pack(C_READ, S_REGB, 8)
my_rom[40] = pack(C_MUL,  S_NONE, 0); my_rom[41] = pack(C_WD,   S_REGC, 12)

# Step 8: t1 = t1 * Sx
my_rom[42] = pack(C_READ, S_REGA, 8); my_rom[43] = pack(C_READ, S_REGB, 11)
my_rom[44] = pack(C_MUL,  S_NONE, 0); my_rom[45] = pack(C_WD,   S_REGC, 8)

# Step 9: t2 = Rx - Rz
my_rom[46] = pack(C_READ, S_REGA, 9); my_rom[47] = pack(C_READ, S_REGB, 12)
my_rom[48] = pack(C_SUB,  S_NONE, 0); my_rom[49] = pack(C_WD,   S_REGC, 10)

# Step 10: Rx = Rx * Rz
my_rom[50] = pack(C_READ, S_REGA, 9); my_rom[51] = pack(C_READ, S_REGB, 12)
my_rom[52] = pack(C_MUL,  S_NONE, 0); my_rom[53] = pack(C_WD,   S_REGC, 9)

# Step 11: Sx2 = A24.x * t2
my_rom[54] = pack(C_READ, S_REGA, 4); my_rom[55] = pack(C_READ, S_REGB, 10)
my_rom[56] = pack(C_MUL,  S_NONE, 0); my_rom[57] = pack(C_WD,   S_REGC, 13)

# Step 12: Sz = t0 - t1
my_rom[58] = pack(C_READ, S_REGA, 7); my_rom[59] = pack(C_READ, S_REGB, 8)
my_rom[60] = pack(C_SUB,  S_NONE, 0); my_rom[61] = pack(C_WD,   S_REGC, 14)

# Step 13: Rz = Rz + Sx2
my_rom[62] = pack(C_READ, S_REGA, 12); my_rom[63] = pack(C_READ, S_REGB, 13)
my_rom[64] = pack(C_ADD,  S_NONE, 0);  my_rom[65] = pack(C_WD,   S_REGC, 12)

# Step 14: Sx = t0 + t1
my_rom[66] = pack(C_READ, S_REGA, 7); my_rom[67] = pack(C_READ, S_REGB, 8)
my_rom[68] = pack(C_ADD,  S_NONE, 0); my_rom[69] = pack(C_WD,   S_REGC, 11)

# Step 15: Rz = Rz * t2
my_rom[70] = pack(C_READ, S_REGA, 12); my_rom[71] = pack(C_READ, S_REGB, 10)
my_rom[72] = pack(C_MUL,  S_NONE, 0);  my_rom[73] = pack(C_WD,   S_REGC, 12)

# Step 16: Sz = Sz^2
my_rom[74] = pack(C_READ, S_REGA, 14); my_rom[75] = pack(C_READ, S_REGB, 14)
my_rom[76] = pack(C_MUL,  S_NONE, 0);  my_rom[77] = pack(C_WD,   S_REGC, 14)

# Step 17: Sx = Sx^2
my_rom[78] = pack(C_READ, S_REGA, 11); my_rom[79] = pack(C_READ, S_REGB, 11)
my_rom[80] = pack(C_MUL,  S_NONE, 0);  my_rom[81] = pack(C_WD,   S_REGC, 11)

# Step 18: Sz = Sz * PQ.x
my_rom[82] = pack(C_READ, S_REGA, 14); my_rom[83] = pack(C_READ, S_REGB, 5)
my_rom[84] = pack(C_MUL,  S_NONE, 0);  my_rom[85] = pack(C_WD,   S_REGC, 14)

# Step 19: Sx = Sx * PQ.z
my_rom[86] = pack(C_READ, S_REGA, 11); my_rom[87] = pack(C_READ, S_REGB, 6)
my_rom[88] = pack(C_MUL,  S_NONE, 0);  my_rom[89] = pack(C_WD,   S_REGC, 11)

# --- Result Write-back ---
my_rom[90] = pack(C_READ, S_REGA, 9);  my_rom[91] = pack(C_WD, S_REGA, 73)
my_rom[92] = pack(C_READ, S_REGA, 12); my_rom[93] = pack(C_WD, S_REGA, 76)
my_rom[94] = pack(C_READ, S_REGA, 11); my_rom[95] = pack(C_WD, S_REGA, 75)
my_rom[96] = pack(C_READ, S_REGA, 14); my_rom[97] = pack(C_WD, S_REGA, 78)

# Termination
my_rom[98] = pack(C_HALT, S_NONE, 0)

# 2. Setup Memory and Run
def hex_to_int(h_str):
    """Helper to convert Verilog hex format to Python int"""
    return int(h_str.replace("255'h", ""), 16)

# ==========================================
# 1. INITIALIZE MEMORY (Mirroring TB execute_ext_cmd)
# ==========================================
initial_mem = {
    # P.x (Addr 64)
    64: Fp2(re=hex_to_int("011dc60f392456de3eb13b9046685257bdd640fb06671ad11c80317fa3b1799d"),
            im=hex_to_int("04b9542316419f828b9d2434e465e150bd9c66b3ad3c2d6d1a3d1fa7bc8960a9")),
    # P.z (Addr 65)
    65: Fp2(re=hex_to_int("04d0ef32815ef6d13b8faa1837f8a88b17fc695a07a0ca6e0822e8f36c031199"),
            im=hex_to_int("035b2d358b8148f6b38a088ca65ed389b74d0fb132e706298fadc1a606cb0fb3")),
    # Q.x (Addr 66)
    66: Fp2(re=hex_to_int("01b8f66b27cd813047229389571aa8766c307511b2b9437a28df6ec4ce4a2bbd"),
            im=hex_to_int("02df309418c267976142ea7d17be31111a2a73ed562b0f79c37459eef50bea63")),
    # Q.z (Addr 67)
    67: Fp2(re=hex_to_int("03ace6f3bacfb3d00b1f9163ce9ff57f43b7a3a69a8dca03580d7b71d8f56413"),
            im=hex_to_int("02586dda8d5288f1142c3fe860e7a113ec1b8ca1f91e1d4c1ff49b7889463e85")),
    # A24.x (Addr 68)
    68: Fp2(re=hex_to_int("048995b1f16287e4e9c349e03602f8ac10f1bc81448aaa9e66b2bc5b50c187fc"),
            im=hex_to_int("04ab134fe5d7b8756dadd6c795a76d79bf3c4c06434308bc89fa6a688fb5d27b")),
    # PQ.x (Addr 69)
    69: Fp2(re=hex_to_int("0189ce9993cd59bf5c941cf0dc98d2c1e2acf72f9e574f7aa0ee89aed453dd32"),
            im=hex_to_int("03a0959447294739614ff3d719db3ad0ddd1dfb23b982ef8daf61a26146d3f31")),
    # PQ.z (Addr 70)
    70: Fp2(re=hex_to_int("0228da676123fdf77656af7229d4beef3eabedcbbaa80dd488bd64072bcfbe01"),
            im=hex_to_int("0298218baf42e12f3838b3268e944239b02b61c4a3d70628ece66fa2fd5166e6"))
}
cu = ControlUnitSim(initial_mem)
engine = RomEngineSim(cu, my_rom)

print("--- EXECUTING PYTHON xDBLADD ---")
engine.run(verbose=False) 

# ==========================================
# 3. VERIFY RESULTS (Mirroring TB $display)
# ==========================================
expected_results = {
    73: ("Rx", "00dfab52916b79c83b49c256ade50de4752da82e7ba79f27516bb40d0618084b", 
               "0357102bfd03a74e4dd5a4a4b1d5d1f8c25573f9907f030d164ed714b2d54781"),
    76: ("Rz", "02b30389e9caa3494baa5c0b1936f06e491786c572c0ef93beed90246fec32f9", 
               "04cf6da21e234794068add1be31f43bc6f98e17597c995a279d7727ff6e96845"),
    75: ("Sx", "044618680d33a4eb843cb979f9ac38dd96e57400d63085f655a5b2bce0910bb5", 
               "03bf282bf26386fcd01cde5ede2a48027acd2093d44639d91a47770d9d280c7d"),
    78: ("Sz", "014be546428342262dd73c8092a90c472e46590fc992f04a3be344d78b003920", 
               "01bd65748a1f31ac8f3ce28f135e4b9e67477f0cfa0de2762408087d81e4dfd8")
}

print("\n--- VERIFYING PYTHON RESULTS ---")
for addr, (name, exp_re, exp_im) in expected_results.items():
    actual = cu.memory.get(addr, Fp2(0, 0))
    exp_re_int = int(exp_re, 16)
    exp_im_int = int(exp_im, 16)
    
    re_pass = "PASS" if actual.re == exp_re_int else "FAIL"
    im_pass = "PASS" if actual.im == exp_im_int else "FAIL"
    
    print(f"{name} (Addr {addr}):")
    print(f"  Real: {hex(actual.re)[2:].zfill(64)} [{re_pass}]")
    print(f"  Imag: {hex(actual.im)[2:].zfill(64)} [{im_pass}]")

print("\nTEST COMPLETE")