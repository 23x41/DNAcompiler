#!/usr/bin/env python3
"""
HolyGATC / GATC46
A genome-sequence compiler and stack VM.
46 opcodes — one for each human chromosome.

MADE THE COMPILER WITH SUPERGROK
"""

from __future__ import annotations

import argparse
import random
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
VERSION = "1.1.0"
BANNER = r"""
  _   _       _        ____    _  _____ ____
 | | | | ___ | |_   _ / ___|  / \|_   _/ ___|
 | |_| |/ _ \| | | | | |  _  / _ \ | || |
 |  _  | (_) | | |_| | |_| |/ ___ \| || |___
 |_| |_|\___/|_|\__, |\____/_/   \_\_| \____|
                |___/
  Genome Assembly & Translation Compiler
  46 opcodes — one for each human chromosome
  MADE THE COMPILER WITH SUPERGROK
"""

GENETIC_CODE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

OPDEF = [
    ("A", "ADD",   0x01, False, "(a b -- a+b)     integer add"),
    ("B", "JMP",   0x02, False, "(addr --)        jump to address"),
    ("C", "CALL",  0x03, False, "(addr --)        call subroutine"),
    ("D", "DUP",   0x04, False, "(a -- a a)       duplicate TOS"),
    ("E", "EQ",    0x05, False, "(a b -- flag)    1 if a == b else 0"),
    ("F", "FETCH", 0x06, False, "(addr -- val)    load memory[addr]"),
    ("G", "GT",    0x07, False, "(a b -- flag)    1 if a > b else 0"),
    ("H", "HALT",  0x08, False, "(--)             stop the machine"),
    ("I", "IN",    0x09, False, "(-- n)           read integer from stdin"),
    ("J", "JZ",    0x0A, False, "(flag addr --)   jump if flag == 0"),
    ("K", "LIT",   0x0B, True,  "(-- n)           push immediate integer"),
    ("L", "LT",    0x0C, False, "(a b -- flag)    1 if a < b else 0"),
    ("M", "MUL",   0x0D, False, "(a b -- a*b)     integer multiply"),
    ("N", "NOT",   0x0E, False, "(a -- flag)      1 if a == 0 else 0"),
    ("O", "OR",    0x0F, False, "(a b -- a|b)     bitwise OR"),
    ("P", "AND",   0x10, False, "(a b -- a&b)     bitwise AND"),
    ("Q", "QUIT",  0x11, False, "(code --)        exit with status code"),
    ("R", "RET",   0x12, False, "(--)             return from CALL"),
    ("S", "SUB",   0x13, False, "(a b -- a-b)     integer subtract"),
    ("T", "STORE", 0x14, False, "(val addr --)    memory[addr] = val"),
    ("U", "OVER",  0x15, False, "(a b -- a b a)   copy NOS to TOS"),
    ("V", "DIV",   0x16, False, "(a b -- a/b)     integer divide"),
    ("W", "PUTC",  0x17, False, "(n --)           write TOS as a character"),
    ("X", "XOR",   0x18, False, "(a b -- a^b)     bitwise XOR"),
    ("Y", "PUTN",  0x19, False, "(n --)           write TOS as a decimal number"),
    ("Z", "SWAP",  0x1A, False, "(a b -- b a)     swap top two"),
    ("_", "NOP",   0x1B, False, "(--)             no operation"),
    ("%", "MOD",   0x1C, False, "(a b -- a%b)     integer modulo"),
    ("~", "NEG",   0x1D, False, "(a -- -a)        negate"),
    ("!", "DROP",  0x1E, False, "(a --)           discard TOS"),
    ("?", "GETC",  0x1F, False, "(-- n)           read one character (ASCII)"),
    ("<", "SHL",   0x20, False, "(a b -- a<<b)    shift left"),
    (">", "SHR",   0x21, False, "(a b -- a>>b)    shift right"),
    ("$", "RAND",  0x22, False, "(n -- r)         random int in [0, n)"),
    ("+", "INC",   0x23, False, "(a -- a+1)       increment"),
    ("-", "DEC",   0x24, False, "(a -- a-1)       decrement"),
    ("*", "JNZ",   0x25, False, "(flag addr --)   jump if flag != 0"),
    ("@", "ROT",   0x26, False, "(a b c -- b c a) rotate top three"),
    ("&", "NIP",   0x27, False, "(a b -- b)       drop second"),
    ("=", "NE",    0x28, False, "(a b -- flag)    1 if a != b else 0"),
    ("/", "GE",    0x29, False, "(a b -- flag)    1 if a >= b else 0"),
    ("^", "LE",    0x2A, False, "(a b -- flag)    1 if a <= b else 0"),
    ("|", "MAX",   0x2B, False, "(a b -- max)     greater of two"),
    (";", "MIN",   0x2C, False, "(a b -- min)     lesser of two"),
    (".", "ABS",   0x2D, False, "(a -- |a|)       absolute value"),
    (",", "DEPTH", 0x2E, False, "(-- n)           push stack depth"),
]

LETTER_TO_OP = {row[0]: row[1] for row in OPDEF}
NAME_TO_BYTE = {row[1]: row[2] for row in OPDEF}
BYTE_TO_NAME = {row[2]: row[1] for row in OPDEF}
BYTE_HAS_IMM = {row[2]: row[3] for row in OPDEF}
LETTER_HAS_IMM = {row[0]: row[3] for row in OPDEF}

ALIASES = {
    "ADD": "ADD", "JMP": "JMP", "CALL": "CALL", "DUP": "DUP", "EQ": "EQ",
    "FETCH": "FETCH", "LOAD": "FETCH", "GT": "GT", "HALT": "HALT", "IN": "IN",
    "JZ": "JZ", "LIT": "LIT", "PUSH": "LIT", "LT": "LT", "MUL": "MUL",
    "NOT": "NOT", "OR": "OR", "AND": "AND", "QUIT": "QUIT", "RET": "RET",
    "SUB": "SUB", "STORE": "STORE", "OVER": "OVER", "DIV": "DIV",
    "PUTC": "PUTC", "EMIT": "PUTC", "XOR": "XOR", "PUTN": "PUTN",
    "PRINT": "PUTN", "SWAP": "SWAP", "NOP": "NOP", "MOD": "MOD",
    "NEG": "NEG", "DROP": "DROP", "GETC": "GETC", "SHL": "SHL",
    "SHR": "SHR", "RAND": "RAND", "INC": "INC", "DEC": "DEC",
    "JNZ": "JNZ", "ROT": "ROT", "NIP": "NIP", "NE": "NE",
    "GE": "GE", "LE": "LE", "MAX": "MAX", "MIN": "MIN",
    "ABS": "ABS", "DEPTH": "DEPTH",
}

MEM_SIZE = 4096
STACK_LIMIT = 65536
CALL_LIMIT = 4096
STEP_LIMIT_DEFAULT = 1_000_000
BITS_OF_BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
BASE_OF_BITS = "ACGT"


class GATCError(Exception):
    pass


@dataclass
class Token:
    kind: str
    value: Union[str, int]
    line: int
    col: int


def tokenize_aa(src: str) -> List[Token]:
    tokens: List[Token] = []
    i = 0
    n = len(src)
    line = 1
    col = 1

    def peek(k: int = 1) -> str:
        j = i + k
        return src[j] if j < n else ""

    while i < n:
        ch = src[i]
        if ch in " \t\r":
            i += 1
            col += 1
            continue
        if ch == "\n":
            i += 1
            line += 1
            col = 1
            continue
        if ch == "#":
            while i < n and src[i] != "\n":
                i += 1
                col += 1
            continue
        if ch == '"':
            start_col = col
            i += 1
            col += 1
            buf: List[str] = []
            while i < n and src[i] != '"':
                if src[i] == "\\" and i + 1 < n:
                    nxt = src[i + 1]
                    buf.append({"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"'}.get(nxt, nxt))
                    i += 2
                    col += 2
                    continue
                if src[i] == "\n":
                    raise GATCError(f"line {line}:{start_col}: unterminated string")
                buf.append(src[i])
                i += 1
                col += 1
            if i >= n or src[i] != '"':
                raise GATCError(f"line {line}:{start_col}: unterminated string")
            i += 1
            col += 1
            tokens.append(Token("STR", "".join(buf), line, start_col))
            continue
        if ch == ":":
            start_col = col
            i += 1
            col += 1
            name = []
            while i < n and (src[i].isalnum() or src[i] == "_"):
                name.append(src[i])
                i += 1
                col += 1
            if not name:
                raise GATCError(f"line {line}:{start_col}: empty label")
            tokens.append(Token("LABEL", "".join(name).upper(), line, start_col))
            continue
        if ch.isdigit() or (ch == "-" and peek(1).isdigit()):
            start_col = col
            sign = 1
            if ch == "-":
                sign = -1
                i += 1
                col += 1
            if i < n and src[i] == "0" and peek(1) in "xX":
                i += 2
                col += 2
                digits = []
                while i < n and src[i] in "0123456789abcdefABCDEF":
                    digits.append(src[i])
                    i += 1
                    col += 1
                if not digits:
                    raise GATCError(f"line {line}:{start_col}: bad hex literal")
                num = sign * int("".join(digits), 16)
            else:
                digits = []
                while i < n and src[i].isdigit():
                    digits.append(src[i])
                    i += 1
                    col += 1
                num = sign * int("".join(digits), 10)
            tokens.append(Token("NUM", num, line, start_col))
            continue
        if ch.isalpha() or ch in LETTER_TO_OP:
            start_col = col
            if ch in LETTER_TO_OP and not (peek(1).isalpha() or peek(1) == "_"):
                tokens.append(Token("OP", ch.upper() if ch.isalpha() else ch, line, start_col))
                i += 1
                col += 1
                continue
            word = []
            while i < n and (src[i].isalnum() or src[i] == "_"):
                word.append(src[i])
                i += 1
                col += 1
            w = "".join(word)
            wu = w.upper()
            if len(w) == 1 and w.upper() in LETTER_TO_OP:
                tokens.append(Token("OP", w.upper(), line, start_col))
            elif wu in ALIASES:
                tokens.append(Token("NAME", ALIASES[wu], line, start_col))
            else:
                tokens.append(Token("IDENT", wu, line, start_col))
            continue
        raise GATCError(f"line {line}:{col}: unexpected character {ch!r}")
    return tokens


@dataclass
class Instr:
    name: str
    imm: Optional[int] = None
    imm_label: Optional[str] = None
    line: int = 0


def assemble_tokens(tokens: List[Token]) -> Tuple[List[Instr], Dict[str, int]]:
    program: List[Instr] = []
    labels: Dict[str, int] = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok.kind == "LABEL":
            if tok.value in labels:
                raise GATCError(f"line {tok.line}: duplicate label :{tok.value}")
            labels[str(tok.value)] = len(program)
            i += 1
            continue
        if tok.kind == "STR":
            for ch in str(tok.value):
                program.append(Instr("LIT", imm=ord(ch), line=tok.line))
                program.append(Instr("PUTC", line=tok.line))
            i += 1
            continue
        if tok.kind == "NUM":
            program.append(Instr("LIT", imm=int(tok.value), line=tok.line))
            i += 1
            continue
        if tok.kind == "IDENT":
            program.append(Instr("LIT", imm_label=str(tok.value), line=tok.line))
            i += 1
            continue
        if tok.kind == "NAME":
            name = str(tok.value)
            imm = None
            imm_label = None
            if name == "LIT":
                if i + 1 >= len(tokens):
                    raise GATCError(f"line {tok.line}: LIT requires an operand")
                nxt = tokens[i + 1]
                if nxt.kind == "NUM":
                    imm = int(nxt.value)
                elif nxt.kind == "IDENT":
                    imm_label = str(nxt.value)
                else:
                    raise GATCError(f"line {tok.line}: LIT requires number or label")
                i += 2
            else:
                i += 1
            program.append(Instr(name, imm=imm, imm_label=imm_label, line=tok.line))
            continue
        if tok.kind == "OP":
            letter = str(tok.value)
            name = LETTER_TO_OP[letter]
            imm = None
            imm_label = None
            if LETTER_HAS_IMM[letter]:
                if i + 1 >= len(tokens):
                    raise GATCError(f"line {tok.line}: {name} requires an operand")
                nxt = tokens[i + 1]
                if nxt.kind == "NUM":
                    imm = int(nxt.value)
                elif nxt.kind == "IDENT":
                    imm_label = str(nxt.value)
                else:
                    raise GATCError(f"line {tok.line}: {name} requires number or label")
                i += 2
            else:
                i += 1
            program.append(Instr(name, imm=imm, imm_label=imm_label, line=tok.line))
            continue
        raise GATCError(f"line {tok.line}: cannot assemble token {tok}")
    return program, labels


def resolve(program: List[Instr], labels: Dict[str, int]) -> List[Instr]:
    out = []
    for ins in program:
        if ins.imm_label is not None:
            if ins.imm_label not in labels:
                raise GATCError(f"line {ins.line}: unknown label :{ins.imm_label}")
            out.append(Instr(ins.name, imm=labels[ins.imm_label], line=ins.line))
        else:
            out.append(ins)
    return out


def encode_bytecode(program: List[Instr]) -> bytes:
    buf = bytearray()
    for ins in program:
        b = NAME_TO_BYTE[ins.name]
        buf.append(b)
        if BYTE_HAS_IMM[b]:
            if ins.imm is None:
                raise GATCError(f"line {ins.line}: {ins.name} missing immediate")
            n = int(ins.imm) & 0xFFFFFFFF
            buf.extend(n.to_bytes(4, "little", signed=False))
    return bytes(buf)


def decode_bytecode(data: bytes) -> List[Instr]:
    program: List[Instr] = []
    i = 0
    while i < len(data):
        b = data[i]
        i += 1
        if b not in BYTE_TO_NAME:
            raise GATCError(f"invalid opcode byte 0x{b:02x} at {i - 1}")
        name = BYTE_TO_NAME[b]
        imm = None
        if BYTE_HAS_IMM[b]:
            if i + 4 > len(data):
                raise GATCError("truncated immediate")
            imm = int.from_bytes(data[i:i + 4], "little", signed=False)
            if imm >= 0x80000000:
                imm -= 0x100000000
            i += 4
        program.append(Instr(name, imm=imm))
    return program


def strip_to_bases(text: str) -> str:
    kept = []
    for line in text.splitlines():
        if line.startswith(">"):
            continue
        kept.append(line)
    blob = "\n".join(kept) if kept else text
    return "".join(ch for ch in blob.upper().replace("U", "T") if ch in "ACGT")


def bytecode_to_dna(bc: bytes) -> str:
    """ATG + 16-bit big-endian length + payload + TAA. Bases are A=00 C=01 G=10 T=11."""
    packet = len(bc).to_bytes(2, "big") + bc
    bases = []
    for byte in packet:
        bases.append(BASE_OF_BITS[(byte >> 6) & 3])
        bases.append(BASE_OF_BITS[(byte >> 4) & 3])
        bases.append(BASE_OF_BITS[(byte >> 2) & 3])
        bases.append(BASE_OF_BITS[byte & 3])
    return "ATG" + "".join(bases) + "TAA"


def dna_to_bytecode(text: str) -> bytes:
    bases = strip_to_bases(text)
    if bases.startswith("ATG"):
        bases = bases[3:]
    if bases.endswith("TAA"):
        bases = bases[:-3]
    if len(bases) % 4 != 0:
        raise GATCError("DNA payload is not a multiple of 4 bases")
    raw = bytearray()
    for i in range(0, len(bases), 4):
        codon = bases[i:i + 4]
        byte = 0
        for ch in codon:
            byte = (byte << 2) | BITS_OF_BASE[ch]
        raw.append(byte)
    if len(raw) < 2:
        raise GATCError("DNA genome missing length")
    ln = int.from_bytes(raw[:2], "big")
    payload = bytes(raw[2:2 + ln])
    if len(payload) != ln:
        raise GATCError(f"DNA length header {ln} but payload is {len(payload)}")
    return payload


def codon_value(codon: str) -> int:
    return (BITS_OF_BASE[codon[0]] << 4) | (BITS_OF_BASE[codon[1]] << 2) | BITS_OF_BASE[codon[2]]


def compile_gene(text: str, require_start: bool = True) -> List[Instr]:
    bases = strip_to_bases(text)
    start = 0
    if require_start:
        idx = bases.find("ATG")
        if idx < 0:
            raise GATCError("no START codon ATG found")
        start = idx
    seq = bases[start:]
    program: List[Instr] = []
    i = 0
    first = True
    while i + 2 < len(seq):
        codon = seq[i:i + 3]
        i += 3
        aa = GENETIC_CODE.get(codon)
        if aa is None:
            raise GATCError(f"unknown codon {codon}")
        if aa == "*":
            program.append(Instr("HALT"))
            break
        if first and aa == "M":
            first = False
            continue
        first = False
        if aa not in LETTER_TO_OP:
            raise GATCError(f"amino acid {aa} from codon {codon} has no opcode")
        name = LETTER_TO_OP[aa]
        if name == "LIT":
            if i + 2 >= len(seq):
                raise GATCError("LIT/Lysine codon missing immediate codon")
            imm_codon = seq[i:i + 3]
            i += 3
            program.append(Instr("LIT", imm=codon_value(imm_codon)))
        else:
            if BYTE_HAS_IMM[NAME_TO_BYTE[name]]:
                raise GATCError(f"{name} needs an immediate the gene encoding cannot supply")
            program.append(Instr(name))
    if not program or program[-1].name != "HALT":
        program.append(Instr("HALT"))
    return program


def compile_source(text: str, filename: str = "") -> List[Instr]:
    lower = filename.lower()
    looks_dna = lower.endswith((".dna", ".gatc", ".fa", ".fasta", ".fna"))
    if not looks_dna:
        bases = strip_to_bases(text)
        letters = [c for c in text.upper() if c.isalpha()]
        if bases.startswith("ATG") and len(bases) >= 12:
            non = [c for c in letters if c not in "ACGTU"]
            if len(non) == 0:
                looks_dna = True
    if looks_dna:
        return decode_bytecode(dna_to_bytecode(text))
    tokens = tokenize_aa(text)
    program, labels = assemble_tokens(tokens)
    return resolve(program, labels)


def write_fasta(dna: str, header: str, width: int = 60) -> str:
    lines = [f">{header}"]
    for i in range(0, len(dna), width):
        lines.append(dna[i:i + width])
    return "\n".join(lines) + "\n"


def translate_coding_sequence(text: str) -> str:
    bases = strip_to_bases(text)
    idx = bases.find("ATG")
    seq = bases[idx:] if idx >= 0 else bases
    letters = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = GENETIC_CODE.get(codon, "?")
        letters.append(aa)
        if aa == "*":
            break
    return "".join(letters)


def disasm(program: List[Instr]) -> str:
    lines = []
    for pc, ins in enumerate(program):
        if ins.imm is not None:
            lines.append(f"{pc:04d}: {ins.name:<6} {ins.imm}")
        else:
            lines.append(f"{pc:04d}: {ins.name}")
    return "\n".join(lines)


class VM:
    def __init__(self, program: List[Instr], stdin=None, stdout=None, mem_size: int = MEM_SIZE):
        self.program = program
        self.pc = 0
        self.stack: List[int] = []
        self.calls: List[int] = []
        self.mem = [0] * mem_size
        self.stdin = stdin or sys.stdin
        self.stdout = stdout or sys.stdout
        self.halted = False
        self.exit_code = 0
        self.steps = 0

    def push(self, n: int) -> None:
        if len(self.stack) >= STACK_LIMIT:
            raise GATCError("stack overflow")
        self.stack.append(int(n))

    def pop(self) -> int:
        if not self.stack:
            raise GATCError(f"stack underflow at pc={self.pc}")
        return self.stack.pop()

    def run(self, max_steps: int = STEP_LIMIT_DEFAULT) -> int:
        n = len(self.program)
        while not self.halted:
            if self.steps >= max_steps:
                raise GATCError(f"step limit {max_steps} exceeded (possible infinite loop)")
            if self.pc < 0 or self.pc >= n:
                raise GATCError(f"pc {self.pc} out of range 0..{n - 1}")
            ins = self.program[self.pc]
            self.steps += 1
            self._exec(ins)
        return self.exit_code

    def _exec(self, ins: Instr) -> None:
        op = ins.name
        if op == "LIT":
            self.push(ins.imm if ins.imm is not None else 0)
            self.pc += 1
        elif op == "ADD":
            b, a = self.pop(), self.pop()
            self.push(a + b)
            self.pc += 1
        elif op == "SUB":
            b, a = self.pop(), self.pop()
            self.push(a - b)
            self.pc += 1
        elif op == "MUL":
            b, a = self.pop(), self.pop()
            self.push(a * b)
            self.pc += 1
        elif op == "DIV":
            b, a = self.pop(), self.pop()
            if b == 0:
                raise GATCError("division by zero")
            self.push(a // b)
            self.pc += 1
        elif op == "MOD":
            b, a = self.pop(), self.pop()
            if b == 0:
                raise GATCError("modulo by zero")
            self.push(a % b)
            self.pc += 1
        elif op == "AND":
            b, a = self.pop(), self.pop()
            self.push(a & b)
            self.pc += 1
        elif op == "OR":
            b, a = self.pop(), self.pop()
            self.push(a | b)
            self.pc += 1
        elif op == "XOR":
            b, a = self.pop(), self.pop()
            self.push(a ^ b)
            self.pc += 1
        elif op == "NOT":
            self.push(1 if self.pop() == 0 else 0)
            self.pc += 1
        elif op == "NEG":
            self.push(-self.pop())
            self.pc += 1
        elif op == "EQ":
            b, a = self.pop(), self.pop()
            self.push(1 if a == b else 0)
            self.pc += 1
        elif op == "LT":
            b, a = self.pop(), self.pop()
            self.push(1 if a < b else 0)
            self.pc += 1
        elif op == "GT":
            b, a = self.pop(), self.pop()
            self.push(1 if a > b else 0)
            self.pc += 1
        elif op == "NE":
            b, a = self.pop(), self.pop()
            self.push(1 if a != b else 0)
            self.pc += 1
        elif op == "GE":
            b, a = self.pop(), self.pop()
            self.push(1 if a >= b else 0)
            self.pc += 1
        elif op == "LE":
            b, a = self.pop(), self.pop()
            self.push(1 if a <= b else 0)
            self.pc += 1
        elif op == "MAX":
            b, a = self.pop(), self.pop()
            self.push(a if a > b else b)
            self.pc += 1
        elif op == "MIN":
            b, a = self.pop(), self.pop()
            self.push(a if a < b else b)
            self.pc += 1
        elif op == "ABS":
            self.push(abs(self.pop()))
            self.pc += 1
        elif op == "DUP":
            if not self.stack:
                raise GATCError("stack underflow on DUP")
            self.push(self.stack[-1])
            self.pc += 1
        elif op == "SWAP":
            b, a = self.pop(), self.pop()
            self.push(b)
            self.push(a)
            self.pc += 1
        elif op == "OVER":
            if len(self.stack) < 2:
                raise GATCError("stack underflow on OVER")
            self.push(self.stack[-2])
            self.pc += 1
        elif op == "DROP":
            self.pop()
            self.pc += 1
        elif op == "NIP":
            b = self.pop()
            self.pop()
            self.push(b)
            self.pc += 1
        elif op == "ROT":
            if len(self.stack) < 3:
                raise GATCError("stack underflow on ROT")
            c = self.stack.pop()
            b = self.stack.pop()
            a = self.stack.pop()
            self.push(b)
            self.push(c)
            self.push(a)
            self.pc += 1
        elif op == "INC":
            self.push(self.pop() + 1)
            self.pc += 1
        elif op == "DEC":
            self.push(self.pop() - 1)
            self.pc += 1
        elif op == "DEPTH":
            self.push(len(self.stack))
            self.pc += 1
        elif op == "SHL":
            b, a = self.pop(), self.pop()
            self.push(a << (b & 31))
            self.pc += 1
        elif op == "SHR":
            b, a = self.pop(), self.pop()
            self.push(a >> (b & 31))
            self.pc += 1
        elif op == "FETCH":
            addr = self.pop()
            if addr < 0 or addr >= len(self.mem):
                raise GATCError(f"fetch out of bounds: {addr}")
            self.push(self.mem[addr])
            self.pc += 1
        elif op == "STORE":
            addr = self.pop()
            val = self.pop()
            if addr < 0 or addr >= len(self.mem):
                raise GATCError(f"store out of bounds: {addr}")
            self.mem[addr] = val
            self.pc += 1
        elif op == "PUTC":
            self.stdout.write(chr(self.pop() & 0xFF))
            self.stdout.flush()
            self.pc += 1
        elif op == "PUTN":
            self.stdout.write(str(self.pop()))
            self.stdout.flush()
            self.pc += 1
        elif op == "IN":
            line = self.stdin.readline()
            if line == "":
                self.push(-1)
            else:
                try:
                    self.push(int(line.strip()))
                except ValueError:
                    self.push(0)
            self.pc += 1
        elif op == "GETC":
            ch = self.stdin.read(1)
            self.push(ord(ch) if ch else -1)
            self.pc += 1
        elif op == "RAND":
            n = self.pop()
            self.push(0 if n <= 0 else random.randrange(n))
            self.pc += 1
        elif op == "JMP":
            self.pc = self.pop()
        elif op == "JZ":
            addr = self.pop()
            flag = self.pop()
            self.pc = addr if flag == 0 else self.pc + 1
        elif op == "JNZ":
            addr = self.pop()
            flag = self.pop()
            self.pc = addr if flag != 0 else self.pc + 1
        elif op == "CALL":
            if len(self.calls) >= CALL_LIMIT:
                raise GATCError("call stack overflow")
            addr = self.pop()
            self.calls.append(self.pc + 1)
            self.pc = addr
        elif op == "RET":
            if not self.calls:
                raise GATCError("return without CALL")
            self.pc = self.calls.pop()
        elif op == "NOP":
            self.pc += 1
        elif op == "HALT":
            self.halted = True
            self.exit_code = 0
        elif op == "QUIT":
            self.exit_code = self.pop()
            self.halted = True
        else:
            raise GATCError(f"unimplemented opcode {op}")


def cmd_opcodes() -> None:
    print(f"HolyGATC opcode table — {len(OPDEF)} opcodes\n")
    print(f"{'Sym':<4} {'Name':<7} {'Byte':<6} {'Imm':<5} Description")
    print("-" * 72)
    for letter, name, byte, imm, desc in OPDEF:
        print(f"{letter:<4} {name:<7} 0x{byte:02X}  {str(imm):<5} {desc}")


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="gatc",
        description="HolyGATC — genome-sequence compiler and VM. 46 opcodes. MADE THE COMPILER WITH SUPERGROK",
    )
    p.add_argument("command", nargs="?", default="help",
                   choices=["run", "compile", "disasm", "translate", "opcodes", "help", "version"])
    p.add_argument("source", nargs="?", help="source file (.aa / .dna / .gatc)")
    p.add_argument("-o", "--output", help="output file")
    p.add_argument("--dna", action="store_true", help="emit packed DNA genome instead of bytecode")
    p.add_argument("--max-steps", type=int, default=STEP_LIMIT_DEFAULT)
    p.add_argument("--no-banner", action="store_true")
    p.add_argument("--gene", action="store_true",
                   help="treat DNA as a coding sequence: ATG..STOP translated to opcodes")
    args = p.parse_intermixed_args(argv)

    if args.command == "help" and not args.source:
        print(BANNER)
        p.print_help()
        print("\nExamples:")
        print("  python GATC46.py run examples/hello.aa")
        print("  python GATC46.py compile examples/hello.aa -o hello.dna --dna")
        print("  python GATC46.py run hello.dna")
        print("  python GATC46.py opcodes")
        return 0
    if args.command == "version":
        print(f"HolyGATC {VERSION} — 46 opcodes — MADE THE COMPILER WITH SUPERGROK")
        return 0
    if args.command == "opcodes":
        cmd_opcodes()
        return 0
    if not args.source:
        print("error: source file required", file=sys.stderr)
        return 2

    try:
        with open(args.source, "r", encoding="utf-8") as f:
            text = f.read()
    except OSError as e:
        print(f"error: cannot read {args.source}: {e}", file=sys.stderr)
        return 1

    try:
        if args.command == "translate":
            out = translate_coding_sequence(text) + "\n"
            if args.output:
                with open(args.output, "w", encoding="utf-8") as f:
                    f.write(out)
            else:
                sys.stdout.write(out)
            return 0

        if args.command == "compile":
            program = compile_gene(text) if args.gene else compile_source(text, args.source)
            bc = encode_bytecode(program)
            if args.dna:
                dna = bytecode_to_dna(bc)
                payload = write_fasta(
                    dna,
                    f"HolyGATC genome from {args.source} | MADE THE COMPILER WITH SUPERGROK",
                )
                if args.output:
                    with open(args.output, "w", encoding="utf-8") as f:
                        f.write(payload)
                else:
                    sys.stdout.write(payload)
            else:
                target = args.output or (args.source + ".gbc")
                with open(target, "wb") as f:
                    f.write(b"GATC")
                    f.write(len(bc).to_bytes(4, "little"))
                    f.write(bc)
                print(f"wrote {target} ({len(bc)} bytes)")
            return 0

        if args.command == "disasm":
            if args.gene:
                program = compile_gene(text)
            elif args.source.lower().endswith(".gbc"):
                raw = open(args.source, "rb").read()
                if raw[:4] == b"GATC":
                    ln = int.from_bytes(raw[4:8], "little")
                    program = decode_bytecode(raw[8:8 + ln])
                else:
                    program = decode_bytecode(raw)
            else:
                program = compile_source(text, args.source)
            sys.stdout.write(disasm(program) + "\n")
            return 0

        if not args.no_banner:
            print(BANNER, file=sys.stderr)
        if args.gene:
            program = compile_gene(text)
        elif args.source.lower().endswith(".gbc"):
            raw = open(args.source, "rb").read()
            if raw[:4] == b"GATC":
                ln = int.from_bytes(raw[4:8], "little")
                program = decode_bytecode(raw[8:8 + ln])
            else:
                program = decode_bytecode(raw)
        else:
            program = compile_source(text, args.source)
        return VM(program).run(max_steps=args.max_steps)
    except GATCError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
