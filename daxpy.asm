section .bss
section .data
section .text
global daxpy_asm

daxpy_asm:
    mov     r10, [rsp+40]       ; r10 = z

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