; ============================================================================
; daxpy.asm
; DAXPY kernel:  Z[i] = A * X[i] + Y[i]
;
; Assemble with NASM targeting Win64 object format, then link with the
; Visual Studio C project:
;
;     nasm -f win64 daxpy.asm -o daxpy.obj
;
; Add daxpy.obj to the Visual Studio project (or run: cl main.c daxpy.obj)
; and declare the prototype in main.c as:
;
;     extern void daxpy_asm(int n, double a, double *x, double *y, double *z);
;
; ----------------------------------------------------------------------------
; Microsoft x64 calling convention (used by MSVC / Visual Studio on Windows):
; arguments are assigned to registers purely by POSITION, and whether that
; position's register is integer or floating-point depends on the argument's
; type:
;
;   arg1 (int n)      -> ECX
;   arg2 (double a)   -> XMM1   (2nd position -> XMM1, not XMM0)
;   arg3 (double *x)  -> R8
;   arg4 (double *y)  -> R9
;   arg5 (double *z)  -> stack, at [rsp+40] on entry
;                        (rsp+0 = return address, rsp+8..rsp+39 = 32-byte
;                         shadow space reserved by the caller, rsp+40 =
;                         first stack-passed argument)
;
; Only scalar SIMD (SSE2) instructions are used: MOVSD, MULSD, ADDSD,
; operating on individual double-precision values in XMM registers.
; ============================================================================

section .text
global daxpy_asm

daxpy_asm:
    ; ---- fetch 5th argument (z pointer) from the stack ----
    mov     r10, [rsp+40]       ; r10 = z

    ; ---- n was passed in ecx; mov ecx,ecx zero-extends into rcx ----
    mov     ecx, ecx            ; rcx = n (zero-extended)
    xor     rax, rax            ; rax = loop index i = 0

    test    rcx, rcx
    jle     .done                ; if n <= 0, nothing to do

.loop:
    movsd   xmm0, [r8 + rax*8]  ; xmm0 = x[i]
    mulsd   xmm0, xmm1          ; xmm0 = a * x[i]
    addsd   xmm0, [r9 + rax*8]  ; xmm0 = a * x[i] + y[i]
    movsd   [r10 + rax*8], xmm0 ; z[i] = xmm0

    inc     rax
    cmp     rax, rcx
    jl      .loop

.done:
    ret