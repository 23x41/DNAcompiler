HolyGATC

Genome Assembly & Translation Compiler

MADE THE COMPILER WITH SUPERGROK

A compiler and stack virtual machine whose instruction set is the alphabet of
life. Every letter A–Z is an opcode. Programs can be written as amino-acid
assembly or stored as a GATC genome.

Inspired by the spirit of TempleOS — one mind, one language, no committee —
and by the standard genetic code (NCBI / codon table). This is recreational
systems programming, not a biology lab protocol.
 | | | | ___ | |_   _ / |  / |   / _|
 | || |/ _ | | | | | |  _  / _ \ | || |
 |  _  | (_) | | || | || |/ ___ | || |_
 || ||_/||_, |__//   __| _|
                |_/

 What you get

- 36 opcodes (A–Z plus 10 extra symbols)
- Amino-acid assembly (`.aa`) with labels, numbers, strings, comments
- Packed DNA genomes (`.dna`) — bytecode encoded as bases `A=00 C=01 G=10 T=11`
- Gene mode — a real coding sequence `ATG … STOP` translates through the
  genetic code into opcodes
- Stack VM: 4096-cell memory, call stack, integer arithmetic, I/O

## Quick start

```bash
python3 gatc.py opcodes
python3 gatc.py run examples/hello.aa
python3 gatc.py run examples/add.aa --no-banner
python3 gatc.py run examples/fib.aa --no-banner
python3 gatc.py compile examples/hello.aa --dna -o hello.dna
python3 gatc.py run hello.dna --no-banner
python3 gatc.py run --gene examples/add.gene.dna --no-banner
python3 gatc.py disasm examples/hello.aa
python3 gatc.py translate examples/add.gene.dna
