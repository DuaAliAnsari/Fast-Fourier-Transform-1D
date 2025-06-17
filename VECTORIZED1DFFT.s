.equ STDOUT, 0xd0580000
.set vectorSize, 1024           # Maximum vector length supported by hardware (number of 32-bit elements)
.set fftSize, 1024              # FFT input size (must be a power of two)

.section .data
realData:       .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 0.960      

imagData:       .float 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
                .float 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8, 1,2,3,4, 5,6,7,8
                .float 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
                .float 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
                .float 0.960     

indexVector:    .space 4096           # 1024 integers (filled by initIndexVector)
bitrev_table:   .space 2048           # 1024 uint16_t (filled by generate_bitrev_table_1024)
fftLength:      .word 1024            # FFT size
twiddleStep:    .float 0.006135923    # 2π/1024
PI:             .float 3.14159265358979323846


.section .text
.global _start
_start:
main:
# Step 1: initialize index vector [0..n-1](:
call initIndexVector

# Step 2: generate bit-reversal table for N=1024
la a0, bitrev_table
call generate_bitrev_table_1024

# Step 2.5: Generate vectorized bit-reversed indices using bit_reversal(:
li a1, 1024                 # N = 1024
la t0, indexVector
vsetvli t1, a1, e32, m1
vle32.v v2, (t0)              # v2 = [0,1,2,...,1023]
call bit_reversal           # v2 now holds bit-reversed indices

# Step 3: Load addresses and FFT size (:
la t0, realData               # Address of real input array
la t1, imagData               # Address of imaginary input array
la t2, fftLength              # Address of FFT size variable
lw t2, 0(t2)                 # Load FFT size value into t2
la t3, twiddleStep            # Address of twiddle factor step
flw fa5, 0(t3)               # Load floating-point twiddle step into fa5

# Step 4: Reorder input arrays according to bit-reversal indices (:
mv a0, t0                    # realData pointer
mv a1, t1                    # imagData pointer
mv a2, t2                    # FFT size
la a3, bitrev_table          # bit-reversal table pointer
call bitReversalReordering

# Step 5: Call your vectorized FFT butterfly routine
mv a0, t0                    # realData pointer
mv a1, t1                    # imagData pointer
mv a2, t2                    # FFT size
call vFFT

# Step 6: Call function to print results
call printResults

# End program infinite loop
j programEnd

#------------------------------------------------------------------------------
# initIndexVector
# Fills indexVector with sequential integers [0..vectorSize-1]
initIndexVector:
la t4, indexVector           # Base address of indexVector
li t5, vectorSize            # Number of elements
li t6, 0                     # Counter i = 0
init_loop:
bge t6, t5, init_done        # exit when i >= vectorSize
sw t6, 0(t4)                 # Store i at indexVector[i]
addi t4, t4, 4               
addi t6, t6, 1               # i++ incrementer
j init_loop
init_done:
    ret

#------------------------------------------------------------------------------
# generate_bitrev_table_1024
# Generates bit-reversal lookup table for N=1024 (10-bit reversal)
# a0 = pointer to bitrev_table uint16_t array for storing in 16 bits 
generate_bitrev_table_1024:
li t0, 1024                  # N = 1024
li t1, 10                    # bits = 10 (since 2^10 = 1024)
li t2, 0                     # i = 0

gen_loop:
bge t2, t0, gen_done        

mv a0, t2                   
mv a1, t1                   
call bit_reverse_10bits     

slli t3, t2, 1              # Calculate offset into bitrev_table: i * 2 bytes (uint16_t)
la t4, bitrev_table          # Load base address of bitrev_table
add t4, t4, t3              # Calculate address of bitrev_table[i]
sh a0, 0(t4)                # Store reversed index (16-bit) at bitrev_table[i]

addi t2, t2, 1             
j gen_loop                  

gen_done:
ret                         


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
# bit_reversal
# Vectorized computation of bit-reversed indices for [0..N-1]
# Input:  v2 = [0,1,2,...,N-1], a1 = N
# Output: v2 = bit-reversed indices
bit_reversal:
addi sp, sp, -20          # Make space on stack for t0-t4
sw t0, 0(sp)
sw t1, 4(sp)
sw t2, 8(sp)
sw t3, 12(sp)
sw t4, 16(sp)

vsetvli t0, a1, e32, m1   # Set vector length for 32-bit elements
vmv.v.i v0, 0            

addi sp, sp, -8          
sw a0, 0(sp)
sw ra, 4(sp)

mv a0, a1                 # a0 = N
call calcLog2Floor        # Compute log2(N), result in a0
addi a0, a0, 1            # log2(N) = floor(log2(N)) + 1 for exact bit count
mv t1, a0                 # t1 = number of bits

lw a0, 0(sp)             
lw ra, 4(sp)             
addi sp, sp, 8

addi t2, zero, 1          

vector_reverse_loop:
bgt t2, t1, vector_reverse_end_loop

sub t3, t1, t2            # t3 = log2(N) - i
li t4, 1
sll t3, t4, t3            # t3 = 1 << log2(N) - i

vand.vx v1, v2, t3        # For each number, mask to check if bit is set
vmsne.vx v3, v1, zero     # v3 = 1 if bit was set, else 0

li t3, 1
sub t4, t2, t3            # t4 = i - 1
sll t4, t3, t4            # t4 = 1 << (i - 1)


vsetvli zero, a1, e32, m1  # Set vector length again to ensure proper config
vmv.v.x v4, t4             # Load the value t4 into v4 
vmerge.vvm v0, v0, v4, v3  # Use vmerge with mask v3 to conditionally add bits

addi t2, t2, 1
j vector_reverse_loop

vector_reverse_end_loop:
vmv.v.v v2, v0            # Copy result (bit-reversed numbers) to v2

lw t0, 0(sp)
lw t1, 4(sp)
lw t2, 8(sp)
lw t3, 12(sp)
lw t4, 16(sp)
addi sp, sp, 20

ret

#------------------------------------------------------------------------------
# bitReversalReordering
# Reorders FFT input arrays (real and imag) according to bit-reversal indices
# Uses RVV vector gather/scatter instructions and precomputed bit-reversal table
#
# a0 = base address of real[]
# a1 = base address of imag[]
# a2 = FFT size N
# a3 = base address of bitrev_table (uint16_t array)
bitReversalReordering:
addi sp, sp, -48            # Allocate 48 bytes on stack for saving registers
sw s0, 0(sp)                # Save callee-saved register s0
sw s1, 4(sp)                
sw s2, 8(sp)                
sw s3, 12(sp)              
sw s4, 16(sp)               
sw ra, 20(sp)              
sw t0, 24(sp)               
sw t1, 28(sp)               
sw t2, 32(sp)               
sw t3, 36(sp)               
sw t4, 40(sp)               
sw t5, 44(sp)              

mv s0, a0                   
mv s1, a1                   
mv s2, a2                  
mv s3, a3                 

vsetvli t0, s2, e32, m1   

la t1, indexVector
vle32.v v0, (t1)             

# Load bit-reversed indices and zero-extend to 32-bit
vle16.v v1, (s3)              # Load 16-bit bit-reversed indices
vzext.vf2 v1, v1              # Zero-extend 16-bit values to 32-bit

vsll.vi v1, v1, 2            


vluxei32.v v2, (s0), v1      # Gather load real data at bit-reversed indices
vluxei32.v v3, (s1), v1      # Gather load imag data at bit-reversed indices


vsll.vi v4, v0, 2          

# Scatter store reordered data
vsuxei32.v v2, (s0), v4      # Scatter store reordered real data
vsuxei32.v v3, (s1), v4      # Scatter store reordered imag data

# Restore registers
lw s0, 0(sp)
lw s1, 4(sp)
lw s2, 8(sp)
lw s3, 12(sp)
lw s4, 16(sp)
lw ra, 20(sp)
lw t0, 24(sp)
lw t1, 28(sp)
lw t2, 32(sp)
lw t3, 36(sp)
lw t4, 40(sp)
lw t5, 44(sp)
addi sp, sp, 48

ret

#------------------------------------------------------------------------------
# calcLog2Floor
# Calculates floor(log2(N)) for a positive integer N.
# This is used to determine the number of fft stages 
#------------------------------------------------------------------------------
calcLog2Floor:
addi sp, sp, -4             
sw t0, 0(sp)                 
mv t0, a0                
li a0, 0                   
logLoop:
beqz t0, logLoopEnd        
srai t0, t0, 1             
addi a0, a0, 1           
j logLoop                  
logLoopEnd:
addi a0, a0, -1            
lw t0, 0(sp)              
addi sp, sp, 4           
ret                      

#------------------------------------------------------------------------------
# vector_sin
# Vectorized sine approximation for input vector in v5
# Output: v5 (overwritten with sine values)
vector_sin:
addi sp,sp, -28           # Allocate stack space
sw a1, 0(sp)
sw t0, 4(sp)
sw t1, 8(sp)
sw t2, 12(sp)
sw t3, 16(sp)
sw t4, 20(sp)
fsw ft0, 24(sp)

vsetvli a1, a0, e32, m1

vmv.v.v v4, v5            # v4 = input vector (angle)
vmv.v.v v6, v5            # v6 = input vector (accumulator)

li t0,-1
li t1,3
li t2,21
li t3,-1

vector_sin_loop:
bgt t1, t2, vector_sin_end_loop

addi t4, t1, -1
mul t4, t4, t1

fcvt.s.w ft0, t4

vfmul.vv v7, v5, v5
vfmul.vv v7, v4, v7
vfdiv.vf v4, v7, ft0

fcvt.s.w ft0, t0

vfmacc.vf v6, ft0, v4     # v6 += t0 * v4

mul t0, t0, t3
addi t1, t1, 2
j vector_sin_loop

vector_sin_end_loop:
vmv.v.v v5, v6

lw a1, 0(sp)
lw t0, 4(sp)
lw t1, 8(sp)
lw t2, 12(sp)
lw t3, 16(sp)
lw t4, 20(sp)
flw ft0, 24(sp)
addi sp, sp, 28
ret

#------------------------------------------------------------------------------
# vector_cos
# Vectorized cosine approximation for input vector in v11
# Output: v11 (overwritten with cosine values)
vector_cos:
addi sp, sp, -28
sw a1, 0(sp)
sw t0, 4(sp)
sw t1, 8(sp)
sw t2, 12(sp)
sw t3, 16(sp)
sw t4, 20(sp)
fsw ft1, 24(sp)

vsetvli a1, a0, e32, m1
vmv.v.i v8, 1
vfcvt.f.x.v v8, v8
vmv.v.i v9, 1
vfcvt.f.x.v v9, v9

li t0, -1
li t1, 2
li t2, 21
li t3, -1

vector_cos_loop:
bgt t1, t2, vector_cos_end_loop
addi t4, t1, -1
mul t4, t4, t1
fcvt.s.w ft1, t4

vfmul.vv v10, v11, v11
vfmul.vv v10, v8, v10
vfdiv.vf v8, v10, ft1
fcvt.s.w ft1, t0
vfmacc.vf v9, ft1, v8     # v9 += t0 * v8

mul t0, t0, t3
addi t1, t1, 2
j vector_cos_loop

vector_cos_end_loop:
vmv.v.v v11, v9

lw a1, 0(sp)
lw t0, 4(sp)
lw t1, 8(sp)
lw t2, 12(sp)
lw t3, 16(sp)
lw t4, 20(sp)
flw ft1, 24(sp)
addi sp, sp, 28
ret

#------------------------------------------------------------------------------
# vFFT
# Vectorized FFT butterfly computation
# a0 = pointer to realData
# a1 = pointer to imagData
# a2 = FFT size (N)
# Assumes input is already bit-reversal reordered!
vFFT:
addi sp, sp, -40
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

mv s0, a0             # s0 = realData pointer
mv s1, a1             # s1 = imagData pointer
mv s2, a2             # s2 = FFT size N

# Calculate log2(N) and store in s3 (number of stages)
mv a0, s2
call calcLog2Floor
addi s3, a0, 1        # s3 = log2(N) (number of stages)

li s4, 1              # s4 = m = 1 (butterfly width)
li s5, 0              # s5 = stage counter

fft_stage_loop:
bge s5, s3, fft_done  # If stage >= log2(N), done

slli s6, s4, 1        # s6 = m2 = m*2 (distance between butterflies)
li t0, 0              # t0 = j (butterfly group start index)

# --- GENERATE TWIDDLE ANGLES ---
vsetvli t1, s4, e32, m1
la t2, indexVector
vle32.v v8, (t2)               # v8 = [0,1,...,m-1]
# Compute angle = -2*pi*k/N for k in [0,m-1]
la t3, twiddleStep
flw ft0, 0(t3)               # ft0 = 2*pi/N
vfmul.vf v9, v8, ft0         # v9 = k * (2*pi/N)
vfneg.v v10, v9              # v10 = -angle (negate each element)

# --- GET COSINE VECTOR ---
vmv.v.v v11, v10             # v11 = angles (for cosine)
mv a0, s4                    # a0 = vector length
call vector_cos              # v11 = cos(-angle)

# --- GET SINE VECTOR ---
vmv.v.v v5, v10              # v5 = angles (for sine)
mv a0, s4
call vector_sin              # v5 = sin(-angle)

# --- BUTTERFLY COMPUTATION ---
fft_butterfly_group_loop:
bge t0, s2, fft_stage_next

# Vector setup: process up to m butterflies at once
vsetvli t1, s4, e32, m1

la t2, indexVector
# Load indices for left and right butterflies
vle32.v v12, (t2)            # v12 = [0,1,...,m-1]
vadd.vx v13, v12, t0         # v13 = left indices = j + [0,1,...,m-1]
vadd.vx v14, v13, s4         # v14 = right indices = left + m

# Convert indices to byte offsets (each float is 4 bytes)
vsll.vi v15, v13, 2          # left_offsets = left * 4
vsll.vi v16, v14, 2          # right_offsets = right * 4

# Load real/imag for left and right using indexed addressing
vluxei32.v v20, (s0), v15    # v20 = real[left]
vluxei32.v v21, (s1), v15    # v21 = imag[left]
vluxei32.v v22, (s0), v16    # v22 = real[right]
vluxei32.v v23, (s1), v16    # v23 = imag[right]

# Twiddle multiplication: (real[right] + j*imag[right]) * (cos - j*sin)
# temp_real = real[right]*cos - imag[right]*sin
vfmul.vv v24, v22, v11       # real[right] * cos
vfnmsac.vv v24, v5, v23      #-imag[right] * sin

# temp_imag = imag[right]*cos + real[right]*sin
vfmul.vv v25, v23, v11       #imag[right] * cos
vfmacc.vv v25, v5, v22       #+real[right] * sin

# Butterfly:
vfadd.vv v26, v20, v24       # real[left]  + temp_real
vfsub.vv v27, v20, v24       # real[left]  - temp_real
vfadd.vv v28, v21, v25       # imag[left]  + temp_imag
vfsub.vv v29, v21, v25       # imag[left]  - temp_imag

# Store results using indexed addressing
vsuxei32.v v26, (s0), v15    # Store real[left]
vsuxei32.v v27, (s0), v16    # Store real[right]
vsuxei32.v v28, (s1), v15    # Store imag[left]
vsuxei32.v v29, (s1), v16    # Store imag[right]

add t0, t0, s6               # Next butterfly group
j fft_butterfly_group_loop

fft_stage_next:
slli s4, s4, 1               # m *= 2
addi s5, s5, 1               # stage++
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
addi sp, sp, 40
ret

#------------------------------------------------------------------------------
# vector_IFFT
# Inverse FFT implementation
# a0 = pointer to realData
# a1 = pointer to imagData
# a2 = FFT size (N)
vector_IFFT:
# Making space on stack for ra,s0-s2
addi sp, sp, -16
sw ra, 0(sp)
sw s0, 4(sp)
sw s1, 8(sp)
sw s2, 12(sp)

# Saving input pointer
mv s0, a0             # s0 = real[] pointer
mv s1, a1             # s1 = imag[] pointer
mv s2, a2             # s2 = N

vsetvli t0, s2, e32, m1
# Conjugate input so that we can use fft as ifft
vle32.v v0, (s1)      # Loading imag[]
vfneg.v v0, v0        # imag[i] = -imag[i]
vse32.v v0, (s1)

# Running fft on conjugated input
mv a0, s0
mv a1, s1
mv a2, s2
call vFFT            # Bit reversal + butterfly stages

# Conjugating output
vle32.v v0, (s1)
vfneg.v v0, v0
vse32.v v0, (s1)

fcvt.s.w ft0, s2     # Converting N to float

# Normalizing the real part by dividing every value by N
vle32.v v1, (s0)
vfdiv.vf v1, v1, ft0
vse32.v v1, (s0)

# Normalizing the imaginary part by dividing every value by N
vle32.v v2, (s1)
vfdiv.vf v2, v2, ft0
vse32.v v2, (s1)

# Restoring from stack
lw ra, 0(sp)
lw s0, 4(sp)
lw s1, 8(sp)
lw s2, 12(sp)
addi sp, sp, 16
ret

#------------------------------------------------------------------------------
# printResults
# Outputs real and imaginary float values to STDOUT
printResults:
addi sp, sp, -24
sw t1, 0(sp)            
sw s1, 4(sp)
sw a0, 8(sp)
sw a1, 12(sp)
fsw ft0, 16(sp)
fsw ft1, 20(sp)

la t1, fftLength          # Loading address of fftLength     
lw s1, 0(t1)              # Loading fftLength value into s1

la a0, realData           # a0 points to real[]
la a1, imagData           # a1 points to imag[]

li t1, 0                  # t1=i i = 0
print_loop:
bge t1, s1, print_end  # If i >= N, end loop

flw ft0, 0(a0)         # Loading real[i] into ft0
flw ft1, 0(a1)         # Loading imag[i] into ft1

li t2, STDOUT          # Base address of STDOUT
fsw ft0, 0(t2)         # Printing real part
fsw ft1, 0(t2)         # Printing imag part

addi a0, a0, 4         # Moving to next real[i] (float=4 bytes)
addi a1, a1, 4         # Moving to next imag[i]
addi t1, t1, 1         # i++
j print_loop

print_end:
lw t1, 0(sp)
lw s1, 4(sp)
lw a0, 8(sp)
lw a1, 12(sp)
flw ft0, 16(sp)
flw ft1, 20(sp)
addi sp, sp, 24
ret
    
programEnd:
j programEnd              # Infinite loop to halt program