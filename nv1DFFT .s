.equ STDOUT, 0xd0580000
.set fftSize, 1024           # FFT input size (must be a power of two)

.section .data
realData:       .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 0.960      # Fill the rest with zeros to make 1024 values

imagData:       .float 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
                .float 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
                .float 0.960      # Fill the rest with zeros to make 1024 values

bitrev_table:   .space 2048           # 1024 uint16_t (filled by generate_bitrev_table_1024)
fftLength:      .word 1024            # FFT size
twiddleStep:    .float 0.006135923    # 2π/1024
PI:             .float 3.14159265358979323846
cosTable:       .space 4096           # Pre-computed cosines for twiddle factors (1024 floats)
sinTable:       .space 4096           # Pre-computed sines for twiddle factors (1024 floats)

.section .text
.global _start
_start:
main:
# Step 1: Generate bit-reversal table for N=1024
la a0, bitrev_table
call generate_bitrev_table_1024

# Step 2: Pre-compute sine and cosine tables for twiddle factors
call generate_twiddle_tables

# Step 3: Load addresses and FFT size
la a0, realData           # Address of real input array
la a1, imagData           # Address of imaginary input array
la a2, fftLength
lw a2, 0(a2)              # a2 = FFT size
la a3, bitrev_table       # Address of bit reversal table

# Step 4: Reorder input arrays according to bit-reversal indices
call bitReversalReordering

# Step 5: Call non-vectorized FFT butterfly routine
la a0, realData
la a1, imagData
la a2, fftLength
lw a2, 0(a2)
call FFT

# Step 6: Call function to print results
call printResults

# End program (infinite loop)
j programEnd

#------------------------------------------------------------------------------
# generate_bitrev_table_1024
# Generates bit-reversal lookup table for N=1024 (10-bit reversal)
# a0 = pointer to bitrev_table (uint16_t array)
generate_bitrev_table_1024:
li t0, 1024                  # N = 1024
li t1, 10                    # bits = 10 (since 2^10 = 1024)
li t2, 0                     # i = 0

gen_loop:
bge t2, t0, gen_done         # Exit loop when i >= N

mv a0, t2                    # Move current index i into a0 (function argument)
mv a1, t1                    # Move number of bits (10) into a1 (function argument)
call bit_reverse_10bits      # Call function to reverse 10 bits of i; result returned in a0

slli t3, t2, 1              # Calculate offset into bitrev_table: i * 2 bytes (uint16_t)
la t4, bitrev_table          # Load base address of bitrev_table
add t4, t4, t3              # Calculate address of bitrev_table[i]
sh a0, 0(t4)                # Store reversed index (16-bit) at bitrev_table[i]

addi t2, t2, 1              # Increment i
j gen_loop                  # Repeat loop

gen_done:
ret                         # Return from function

#------------------------------------------------------------------------------
# bit_reverse_10bits
# Reverses the lowest 10 bits of input in a0
# Returns reversed bits in a0
bit_reverse_10bits:
li t0, 0                    # rev = 0
li t1, 0                    # bit counter

bitrev_loop_10:
slti t3, t1, 10             # Set t3 = 1 if t1 < 10, else 0
beqz t3, bitrev_done_10     # Stop after 10 bits (fixed)
slli t0, t0, 1              # rev <<= 1
andi t2, a0, 1              # Extract LSB of a0
or t0, t0, t2               # rev |= LSB
srli a0, a0, 1              # a0 >>= 1
addi t1, t1, 1              # bit counter++
j bitrev_loop_10

bitrev_done_10:
mv a0, t0                   # Return reversed bits in a0
ret

#------------------------------------------------------------------------------
# calcLog2Floor
# Calculates floor(log2(N)) for a positive integer N.
# This is used to determine the number of FFT stages (levels).
#------------------------------------------------------------------------------
calcLog2Floor:
addi sp, sp, -4              # Allocate stack space for t0
sw t0, 0(sp)                 # Save t0 on stack
mv t0, a0                   # Copy input N to t0
li a0, 0                    # Initialize counter to zero
logLoop:
beqz t0, logLoopEnd         # Exit if t0 == 0
srai t0, t0, 1             # Divide t0 by 2 (arithmetic right shift)
addi a0, a0, 1             # Increment division count
j logLoop                  # Repeat loop
logLoopEnd:
addi a0, a0, -1            # Adjust count to floor(log2(N))
lw t0, 0(sp)               # Restore t0
addi sp, sp, 4             # Free stack space
ret                      # Return with result in a0

#------------------------------------------------------------------------------
# generate_twiddle_tables
# Generates sine and cosine tables for all twiddle factors needed in FFT
generate_twiddle_tables:
addi sp, sp, -28
sw ra, 0(sp)
sw s0, 4(sp)
sw s1, 8(sp)
sw s2, 12(sp)
sw s3, 16(sp)
fsw fs0, 20(sp)
fsw fs1, 24(sp)

la s0, cosTable            # s0 = cosTable address
la s1, sinTable            # s1 = sinTable address
la t0, fftLength
lw s2, 0(t0)               # s2 = FFT size (N)
la t0, twiddleStep
flw fs0, 0(t0)            # fs0 = 2π/N

li s3, 0                  # s3 = k (index counter)

twiddle_loop:
bge s3, s2, twiddle_done  # If k >= N, we're done

# Compute angle = -2π * k / N
fcvt.s.w fs1, s3         # Convert k to float
fmul.s fs1, fs1, fs0     # fs1 = k * (2π/N)
fneg.s fs1, fs1          # fs1 = -angle

# Calculate cos(-angle) and sin(-angle)
fmv.s fa0, fs1           # fa0 = angle argument
call scalar_cos          # Returns cos(angle) in fa0
slli t0, s3, 2           # t0 = k * 4 (byte offset for float)
add t1, s0, t0           # t1 = address of cosTable[k]
fsw fa0, 0(t1)           # Store cosine value

fmv.s fa0, fs1           # fa0 = angle argument
call scalar_sin          # Returns sin(angle) in fa0
slli t0, s3, 2           # t0 = k * 4 (byte offset for float)
add t1, s1, t0           # t1 = address of sinTable[k]
fsw fa0, 0(t1)           # Store sine value

addi s3, s3, 1           # k++
j twiddle_loop

twiddle_done:
lw ra, 0(sp)
lw s0, 4(sp)
lw s1, 8(sp)
lw s2, 12(sp)
lw s3, 16(sp)
flw fs0, 20(sp)
flw fs1, 24(sp)
addi sp, sp, 28
ret

#------------------------------------------------------------------------------
# scalar_sin
# Sine approximation using Taylor series
# Input: fa0 = angle in radians
# Output: fa0 = sin(angle)
scalar_sin:
addi sp, sp, -16
fsw fs0, 0(sp)
fsw fs1, 4(sp)
fsw fs2, 8(sp)
fsw fs3, 12(sp)

fmv.s fs0, fa0     # fs0 = x (angle)
fmv.s fs1, fa0     # fs1 = x³/3! + x⁵/5! + ...
fmv.s fs2, fa0     # fs2 = accumulator for result

# Use Taylor series: sin(x) = x - x³/3! + x⁵/5! - ...

# Calculate x³/3!
fmul.s fs3, fs0, fs0   # fs3 = x²
fmul.s fs3, fs3, fs0   # fs3 = x³
li t0, 6               # t0 = 3! = 6
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x³/3!
fneg.s fs3, fs3        # fs3 = -x³/3!
fadd.s fs2, fs2, fs3   # fs2 += -x³/3!

# Calculate x⁵/5!
fmul.s fs3, fs3, fs0   # Reuse previous term, multiply by x²
fmul.s fs3, fs3, fs0
li t0, -20             # t0 = -5!/3! = -20
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x⁵/5!
fadd.s fs2, fs2, fs3   # fs2 += x⁵/5!

# Calculate x⁷/7!
fmul.s fs3, fs3, fs0   # Reuse previous term, multiply by x²
fmul.s fs3, fs3, fs0
li t0, -42             # t0 = -7!/5! = -42
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x⁷/7!
fadd.s fs2, fs2, fs3   # fs2 += x⁷/7!

# Calculate x⁹/9!
fmul.s fs3, fs3, fs0   # Reuse previous term, multiply by x²
fmul.s fs3, fs3, fs0
li t0, -72             # t0 = -9!/7! = -72
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x⁹/9!
fadd.s fs2, fs2, fs3   # fs2 += x⁹/9!

fmv.s fa0, fs2         # Return result in fa0

flw fs0, 0(sp)
flw fs1, 4(sp)
flw fs2, 8(sp)
flw fs3, 12(sp)
addi sp, sp, 16
ret

#------------------------------------------------------------------------------
# scalar_cos
# Cosine approximation using Taylor series
# Input: fa0 = angle in radians
# Output: fa0 = cos(angle)
scalar_cos:
addi sp, sp, -16
fsw fs0, 0(sp)
fsw fs1, 4(sp)
fsw fs2, 8(sp)
fsw fs3, 12(sp)

fmv.s fs0, fa0     # fs0 = x (angle)
li t0, 1
fcvt.s.w fs1, t0   # fs1 = 1.0
fmv.s fs2, fs1     # fs2 = accumulator for result

# Use Taylor series: cos(x) = 1 - x²/2! + x⁴/4! - ...

# Calculate -x²/2!
fmul.s fs3, fs0, fs0   # fs3 = x²
li t0, 2               # t0 = 2! = 2
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x²/2!
fneg.s fs3, fs3        # fs3 = -x²/2!
fadd.s fs2, fs2, fs3   # fs2 += -x²/2!

# Calculate x⁴/4!
fmul.s fs3, fs3, fs0   # Reuse previous term, multiply by x²
fmul.s fs3, fs3, fs0
li t0, -12             # t0 = -4!/2! = -12
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x⁴/4!
fadd.s fs2, fs2, fs3   # fs2 += x⁴/4!

# Calculate -x⁶/6!
fmul.s fs3, fs3, fs0   # Reuse previous term, multiply by x²
fmul.s fs3, fs3, fs0
li t0, -30             # t0 = -6!/4! = -30
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = -x⁶/6!
fadd.s fs2, fs2, fs3   # fs2 += -x⁶/6!

# Calculate x⁸/8!
fmul.s fs3, fs3, fs0   # Reuse previous term, multiply by x²
fmul.s fs3, fs3, fs0
li t0, -56             # t0 = -8!/6! = -56
fcvt.s.w fa0, t0
fdiv.s fs3, fs3, fa0   # fs3 = x⁸/8!
fadd.s fs2, fs2, fs3   # fs2 += x⁸/8!

fmv.s fa0, fs2         # Return result in fa0

flw fs0, 0(sp)
flw fs1, 4(sp)
flw fs2, 8(sp)
flw fs3, 12(sp)
addi sp, sp, 16
ret

#------------------------------------------------------------------------------
# bitReversalReordering
# Reorders FFT input arrays (real and imag) according to bit-reversal indices
#
# a0 = base address of real[]
# a1 = base address of imag[]
# a2 = FFT size N
# a3 = base address of bitrev_table (uint16_t array)
bitReversalReordering:
addi sp, sp, -36           # Allocate stack space for registers
sw ra, 0(sp)
sw s0, 4(sp)
sw s1, 8(sp)
sw s2, 12(sp)
sw s3, 16(sp)
sw s4, 20(sp)
sw s5, 24(sp)
fsw fs0, 28(sp)
fsw fs1, 32(sp)

mv s0, a0                 # s0 = real array base
mv s1, a1                 # s1 = imag array base
mv s2, a2                 # s2 = N (FFT size)
mv s3, a3                 # s3 = bit reversal table

li s4, 0                  # s4 = i (loop counter)

bitrev_loop:
bge s4, s2, bitrev_done   # If i >= N, we're done

# Get bit-reversed index
slli t0, s4, 1            # t0 = i * 2 (byte offset in bitrev_table for uint16)
add t0, s3, t0            # t0 = address of bitrev_table[i]
lhu t1, 0(t0)             # t1 = j = bit-reversed index

# Only swap if j > i (to avoid swapping twice)
ble t1, s4, bitrev_skip

# Load real[i] and imag[i]
slli t0, s4, 2            # t0 = i * 4 (byte offset for float)
add t0, s0, t0            # t0 = address of real[i]
flw fs0, 0(t0)            # fs0 = real[i]
slli t2, s4, 2
add t2, s1, t2            # t2 = address of imag[i]
flw fs1, 0(t2)            # fs1 = imag[i]

# Load real[j] and imag[j]
slli t3, t1, 2            # t3 = j * 4 (byte offset for float)
add t3, s0, t3            # t3 = address of real[j]
flw ft0, 0(t3)            # ft0 = real[j]
slli t4, t1, 2
add t4, s1, t4            # t4 = address of imag[j]
flw ft1, 0(t4)            # ft1 = imag[j]

# Swap real[i] with real[j] and imag[i] with imag[j]
fsw ft0, 0(t0)            # real[i] = real[j]
fsw ft1, 0(t2)            # imag[i] = imag[j]
fsw fs0, 0(t3)            # real[j] = real[i]
fsw fs1, 0(t4)            # imag[j] = imag[i]

bitrev_skip:
addi s4, s4, 1            # i++
j bitrev_loop

bitrev_done:
lw ra, 0(sp)
lw s0, 4(sp)
lw s1, 8(sp)
lw s2, 12(sp)
lw s3, 16(sp)
lw s4, 20(sp)
lw s5, 24(sp)
flw fs0, 28(sp)
flw fs1, 32(sp)
addi sp, sp, 36
ret

#------------------------------------------------------------------------------
# FFT
# Non-vectorized FFT butterfly computation
# a0 = pointer to realData
# a1 = pointer to imagData
# a2 = FFT size (N)
# Assumes input is already bit-reversal reordered!
FFT:
addi sp, sp, -48
sw ra, 0(sp)
sw s0, 4(sp)
sw s1, 8(sp)
sw s2, 12(sp)
sw s3, 16(sp)
sw s4, 20(sp)
sw s5, 24(sp)
sw s6, 28(sp)
sw s7, 32(sp)
sw s8, 36(sp)
sw s9, 40(sp)
sw s10, 44(sp)

mv s0, a0             # s0 = realData pointer
mv s1, a1             # s1 = imagData pointer
mv s2, a2             # s2 = FFT size N

# Calculate log2(N) and store in s3 (number of stages)
mv a0, s2
call calcLog2Floor
addi s3, a0, 1        # s3 = log2(N) (number of stages)

li s4, 1              # s4 = m = 1 (butterfly width)
li s5, 0              # s5 = stage counter

la s9, cosTable       # s9 = pointer to cosine table
la s10, sinTable      # s10 = pointer to sine table

fft_stage_loop:
bge s5, s3, fft_done  # If stage >= log2(N), done

slli s6, s4, 1        # s6 = m2 = m*2 (distance between butterflies)
li s7, 0              # s7 = j (butterfly group start index)

fft_butterfly_group_loop:
bge s7, s2, fft_stage_next

li s8, 0             # s8 = k (offset within group)

fft_butterfly_loop:
bge s8, s4, fft_butterfly_group_next

# Calculate indices for the butterfly operation
add t0, s7, s8       # t0 = idx_left = j + k
add t1, t0, s4       # t1 = idx_right = j + k + m

# Convert indices to byte offsets (each float is 4 bytes)
slli t2, t0, 2       # t2 = left_offset = left * 4
slli t3, t1, 2       # t3 = right_offset = right * 4

# Load real and imaginary parts from memory
add t4, s0, t2       # t4 = address of real[left]
flw ft0, 0(t4)       # ft0 = real[left]
add t5, s1, t2       # t5 = address of imag[left]
flw ft1, 0(t5)       # ft1 = imag[left]
add t6, s0, t3       # t6 = address of real[right]
flw ft2, 0(t6)       # ft2 = real[right]
add t0, s1, t3       # t0 = address of imag[right]
flw ft3, 0(t0)       # ft3 = imag[right]

# Get twiddle factors (cos, sin) for this butterfly
slli t1, s8, 2       # t1 = k * 4 (byte offset in twiddle tables)
mv t2, s2            # t2 = N
div t1, t1, s6       # t1 = (k * N) / (2 * m) (appropriate twiddle index)
add t3, s9, t1       # t3 = address of cosTable[k*N/(2*m)]
flw ft4, 0(t3)       # ft4 = cos
add t3, s10, t1      # t3 = address of sinTable[k*N/(2*m)]
flw ft5, 0(t3)       # ft5 = sin

# Twiddle multiplication: (real[right] + j*imag[right]) * (cos - j*sin)
# temp_real = real[right]*cos - imag[right]*sin
fmul.s ft6, ft2, ft4  # real[right] * cos
fmul.s ft7, ft3, ft5  # imag[right] * sin
fsub.s ft6, ft6, ft7  # temp_real = real[right]*cos - imag[right]*sin

# temp_imag = imag[right]*cos + real[right]*sin
fmul.s ft7, ft3, ft4  # imag[right] * cos
fmul.s ft8, ft2, ft5  # real[right] * sin
fadd.s ft7, ft7, ft8  # temp_imag = imag[right]*cos + real[right]*sin

# Butterfly operation (in-place)
# real[left] = real[left] + temp_real
# real[right] = real[left] - temp_real
# imag[left] = imag[left] + temp_imag
# imag[right] = imag[left] - temp_imag

fadd.s ft8, ft0, ft6  # new_real_left = real[left] + temp_real
fsub.s ft9, ft0, ft6  # new_real_right = real[left] - temp_real
fadd.s ft10, ft1, ft7 # new_imag_left = imag[left] + temp_imag
fsub.s ft11, ft1, ft7 # new_imag_right = imag[left] - temp_imag

# Store results back to memory
fsw ft8, 0(t4)        # real[left] = new_real_left
fsw ft9, 0(t6)        # real[right] = new_real_right
fsw ft10, 0(t5)       # imag[left] = new_imag_left
fsw ft11, 0(t0)       # imag[right] = new_imag_right

addi s8, s8, 1        # k++
j fft_butterfly_loop

fft_butterfly_group_next:
add s7, s7, s6        # j += m2
j fft_butterfly_group_loop

fft_stage_next:
slli s4, s4, 1        # m *= 2
addi s5, s5, 1        # stage++
j fft_stage_loop

fft_done:
lw ra, 0(sp)
lw s0, 4(sp)
lw s1, 8(sp)
lw s2, 12(sp)
lw s3, 16(sp)
lw s4, 20(sp)
lw s5, 24(sp)
lw s6, 28(sp)
lw s7, 32(sp)
lw s8, 36(sp)
lw s9, 40(sp)
lw s10, 44(sp)
addi sp, sp, 48
ret

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# IFFT
# Inverse FFT implementation
# a0 = pointer to realData
# a1 = pointer to imagData
# a2 = FFT size (N)
IFFT:
addi sp, sp, -36
sw ra, 0(sp)
sw s0, 4(sp)
sw s1, 8(sp)
sw s2, 12(sp)
sw t0, 16(sp)
sw t1, 20(sp)
fsw fs0, 24(sp)
fsw fs1, 28(sp)
fsw fs2, 32(sp)

# Save input pointers
mv s0, a0             # s0 = real[] pointer
mv s1, a1             # s1 = imag[] pointer
mv s2, a2             # s2 = N

# Step 1: Conjugate input data (negate imaginary part)
li t0, 0              # Initialize loop counter

ifft_conj_loop1:
bge t0, s2, ifft_conj_done1

# Load imaginary value
slli t1, t0, 2        # t1 = i * 4 (byte offset)
add t1, s1, t1        # t1 = address of imag[i]
flw fs0, 0(t1)        # fs0 = imag[i]

# Negate and store back
fneg.s fs0, fs0       # fs0 = -imag[i]
fsw fs0, 0(t1)        # imag[i] = -imag[i]

addi t0, t0, 1        # i++
j ifft_conj_loop1

ifft_conj_done1:
# Step 2: Perform forward FFT on conjugated data
mv a0, s0             # a0 = realData pointer
mv a1, s1             # a1 = imagData pointer
mv a2, s2             # a2 = FFT size N
call FFT              # Use the existing FFT implementation

# Step 3: Conjugate the result and scale by 1/N
li t0, 0              # Reset counter
fcvt.s.w fs2, s2      # Convert N to float (fs2 = N)
li t1, 1
fcvt.s.w fs1, t1      # fs1 = 1.0
fdiv.s fs1, fs1, fs2  # fs1 = 1/N (scaling factor)

ifft_conj_loop2:
bge t0, s2, ifft_conj_done2

# Get byte offset for this element
slli t1, t0, 2        # t1 = i * 4 (byte offset)

# Load, scale, and store real part
add t1, s0, t1        # t1 = address of real[i]
flw fs0, 0(t1)        # fs0 = real[i]
fmul.s fs0, fs0, fs1  # fs0 = real[i] * (1/N)
fsw fs0, 0(t1)        # real[i] = real[i] * (1/N)

# Load, negate, scale, and store imaginary part
slli t1, t0, 2        # t1 = i * 4 (byte offset)
add t1, s1, t1        # t1 = address of imag[i]
flw fs0, 0(t1)        # fs0 = imag[i]
fneg.s fs0, fs0       # fs0 = -imag[i] (conjugate again)
fmul.s fs0, fs0, fs1  # fs0 = -imag[i] * (1/N)
fsw fs0, 0(t1)        # imag[i] = -imag[i] * (1/N)

addi t0, t0, 1        # i++
j ifft_conj_loop2

ifft_conj_done2:
# Restore saved registers and return
lw ra, 0(sp)
lw s0, 4(sp)
lw s1, 8(sp)
lw s2, 12(sp)
lw t0, 16(sp)
lw t1, 20(sp)
flw fs0, 24(sp)
flw fs1, 28(sp)
flw fs2, 32(sp)
addi sp, sp, 36
ret
#------------------------------------------------------------------------------
# printResults
# Outputs real and imaginary float values to STDOUT
#------------------------------------------------------------------------------
printResults:
addi sp, sp, -24          # Allocate stack space for 6 registers
sw ra, 0(sp)              # Save return address
sw s1, 4(sp)              # Save s1 (callee-saved)
sw t1, 8(sp)              # Save t1 (temporary)
sw t2, 12(sp)             # Save t2 (temporary)
fsw fs0, 16(sp)           # Save fs0 (float temporary)
fsw fs1, 20(sp)           # Save fs1 (float temporary)

la t0, fftLength          # Load address of fftLength
lw s1, 0(t0)              # s1 = fftLength (N)

la t1, realData           # t1 = pointer to real data array
la t2, imagData           # t2 = pointer to imaginary data array

li t0, 0                  # t0 = loop counter (i = 0)

print_loop:
bge t0, s1, print_done    # Exit loop if i >= N

# Calculate byte offsets (i * 4) for current elements
slli t3, t0, 2            # t3 = i * 4 (byte offset)

# Load current values
add t4, t1, t3            # t4 = address of real[i]
flw fs0, 0(t4)            # fs0 = real[i]

add t4, t2, t3            # t4 = address of imag[i]
flw fs1, 0(t4)            # fs1 = imag[i]

# Output to STDOUT
li t4, STDOUT             # t4 = STDOUT memory-mapped address
fsw fs0, 0(t4)            # Write real part to STDOUT
fsw fs1, 0(t4)            # Write imaginary part to STDOUT

addi t0, t0, 1            # Increment loop counter (i++)
j print_loop              # Repeat loop

print_done:
# Restore saved registers
lw ra, 0(sp)              # Restore return address
lw s1, 4(sp)              # Restore s1
lw t1, 8(sp)              # Restore t1
lw t2, 12(sp)             # Restore t2
flw fs0, 16(sp)           # Restore fs0
flw fs1, 20(sp)           # Restore fs1
addi sp, sp, 24           # Free stack space

ret                       # Return to caller

#------------------------------------------------------------------------------
# Program end handler - infinite loop to keep program from running past end
#------------------------------------------------------------------------------
programEnd:
j programEnd              # Infinite loop to halt program