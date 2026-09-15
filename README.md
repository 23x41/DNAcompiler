1. Amino-acid assembly (.aa)
2. Each letter is an opcode.
3.  K (Lysine / LIT) pushes a number.
M is multiply.
H is halt.
Strings emit themselves.



# DNAcompiler

# 40 + 2


K 40

K 2

A

Y

K 10 W

H


2. Packed DNA (.dna)compile --dna writes a FASTA file:

ATG + 16-bit length + payload + TAA

3. Gene mode (--gene)A coding sequence is read like a ribosome:Find ATG (START / Methionine). The start codon is skipped.
Each following codon becomes its amino-acid opcode.
After Lysine (AAA / AAG = K / LIT) the next codon is a number 0–63
in base-4 (A=0 C=1 G=2 T=3).
TAA / TAG / TGA become HALT.

Example gene that prints 42:

ATG AAA GGA AAA AAG GCT TAT TAA
     K   40   K    2  ADD PUTN STOP


40 is codon GGA because G=2,G=2,A=0 → 2*16 + 2*4 + 0 = 40.

LanguageComments start with #. Labels are :NAME. Bare names push that address.
Jumps take the address from the stack.Want
Pattern
Push number
K 7 or just 7
Add / sub / mul
A S M
Print number
Y
Print char
W
Jump
K LABEL then B
Jump if zero
push flag, then K LABEL then J
Call / return
K SUB C … :SUB … R
Memory
T store (val addr --) / F fetch (addr -- val)

JZ stack effect is (flag addr --). It jumps when flag == 0.Put subroutines after a jump to main, or execution will fall into them.

K MAIN
B
:DOUBLE
  D A
  R
:MAIN
  K 21
  K DOUBLE
  C
  Y
  H
