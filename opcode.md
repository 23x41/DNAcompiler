Full table: OPCODES.mdAmino acid → opcodeLetter
Amino acid
Opcode
Letter
Amino acid
Opcode
A
Alanine
ADD
N
Asparagine
NOT
C
Cysteine
CALL
P
Proline
AND
D
Aspartic acid
DUP
Q
Glutamine
QUIT
E
Glutamic acid
EQ
R
Arginine
RET
F
Phenylalanine
FETCH
S
Serine
SUB
G
Glycine
GT
T
Threonine
STORE
H
Histidine
HALT
V
Valine
DIV
I
Isoleucine
IN
W
Tryptophan
PUTC
K
Lysine
LIT
Y
Tyrosine
PUTN
L
Leucine
LT

M
Methionine
MUL

Synthetic letters (not in the natural 20): B JMP, J JZ, O OR, U OVER, X XOR, Z SWAP.START codon ATG is program entry in gene mode. STOP codons TAA TAG TGA halt.



CommandsCommand
Meaning
run FILE
compile if needed and execute
compile FILE -o OUT
write .gbc bytecode
compile FILE --dna -o OUT.dna
write packed FASTA genome
disasm FILE
print instruction listing
translate FILE.dna
DNA → amino-acid letters
run --gene FILE.dna
execute a coding sequence
opcodes
print the 36-opcode table

