#!/usr/bin/env python3
"""HolyGATC / GATC52 version 3.0.0

52 opcodes, 64 codons, 12 third-base synonyms.
Packed DNA: ATG + 16-bit big-endian length + payload + TAA.
A biological FASTA is not packed. run on .dna/.fa/.fasta uses codon
mode unless the file is a real packed genome (ATG, TAA, length matches).
Length mismatch, unknown bytes, return without CALL, and stack underflow
are reported and ignored. The machine runs anyway.
"""
from __future__ import annotations
import argparse, random, sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

VERSION = "3.0.0"
BANNER = r"""
  ____    _  _____ ____ ____ ____
 / ___|  / \|_   _/ ___| ___|___ \
| |  _  / _ \ | || |   |___ \ __) |
| |_| |/ ___ \| || |___ ___) / __/
 \____/_/   \_\_| \____|____/_____|
  Genome Assembly & Translation Compiler
  52 opcodes. 64 codons. Twelve synonyms.
  Version 3. Every codon is an opcode.
"""
OPDEF = [
    ("A","ADD",0x01,False,"(a b -- a+b)","integer add"),
    ("B","JMP",0x02,False,"(addr --)","jump"),
    ("C","CALL",0x03,False,"(addr --)","call"),
    ("D","DUP",0x04,False,"(a -- a a)","duplicate"),
    ("E","EQ",0x05,False,"(a b -- flag)","equal"),
    ("F","FETCH",0x06,False,"(addr -- val)","load"),
    ("G","GT",0x07,False,"(a b -- flag)","greater"),
    ("H","HALT",0x08,False,"(--)","stop"),
    ("I","IN",0x09,False,"(-- n)","read integer"),
    ("J","JZ",0x0A,False,"(flag addr --)","jump if 0"),
    ("K","LIT",0x0B,True,"(-- n)","push immediate"),
    ("L","LT",0x0C,False,"(a b -- flag)","less"),
    ("M","MUL",0x0D,False,"(a b -- a*b)","multiply"),
    ("N","NOT",0x0E,False,"(a -- flag)","1 if zero"),
    ("O","OR",0x0F,False,"(a b -- a|b)","bitwise or"),
    ("P","AND",0x10,False,"(a b -- a&b)","bitwise and"),
    ("Q","QUIT",0x11,False,"(code --)","exit with status"),
    ("R","RET",0x12,False,"(--)","return"),
    ("S","SUB",0x13,False,"(a b -- a-b)","subtract"),
    ("T","STORE",0x14,False,"(val addr --)","store"),
    ("U","OVER",0x15,False,"(a b -- a b a)","copy second"),
    ("V","DIV",0x16,False,"(a b -- a/b)","divide"),
    ("W","PUTC",0x17,False,"(n --)","write char"),
    ("X","XOR",0x18,False,"(a b -- a^b)","xor"),
    ("Y","PUTN",0x19,False,"(n --)","write number"),
    ("Z","SWAP",0x1A,False,"(a b -- b a)","swap"),
    ("_","NOP",0x1B,False,"(--)","nop"),
    ("%","MOD",0x1C,False,"(a b -- a%b)","modulo"),
    ("~","NEG",0x1D,False,"(a -- -a)","negate"),
    ("!","DROP",0x1E,False,"(a --)","drop"),
    ("?","GETC",0x1F,False,"(-- n)","read char"),
    ("<","SHL",0x20,False,"(a b -- a<<b)","shift left"),
    (">","SHR",0x21,False,"(a b -- a>>b)","shift right"),
    ("$","RAND",0x22,False,"(n -- r)","random"),
    ("+","INC",0x23,False,"(a -- a+1)","increment"),
    ("-","DEC",0x24,False,"(a -- a-1)","decrement"),
    ("*","JNZ",0x25,False,"(flag addr --)","jump if not 0"),
    ("@","ROT",0x26,False,"(a b c -- b c a)","rotate"),
    ("&","NIP",0x27,False,"(a b -- b)","nip"),
    ("=","NE",0x28,False,"(a b -- flag)","not equal"),
    ("/","GE",0x29,False,"(a b -- flag)","greater or equal"),
    ("^","LE",0x2A,False,"(a b -- flag)","less or equal"),
    ("|","MAX",0x2B,False,"(a b -- max)","max"),
    (";","MIN",0x2C,False,"(a b -- min)","min"),
    (".","ABS",0x2D,False,"(a -- |a|)","absolute"),
    (",","DEPTH",0x2E,False,"(-- n)","depth"),
    ("`","TUCK",0x2F,False,"(a b -- b a b)","tuck"),
    ("'","PICK",0x30,False,"(n -- x)","pick, 0 is top"),
    ("[","DDUP",0x31,False,"(a b -- a b a b)","2dup"),
    ("]","DDROP",0x32,False,"(a b --)","2drop"),
    ("{","TOR",0x33,False,"(a --)","to return stack"),
    ("}","RFROM",0x34,False,"(-- a)","from return stack"),
]
LETTER_TO_OP = {s: n for s, n, *_ in OPDEF}
NAME_TO_BYTE = {n: b for _, n, b, *_ in OPDEF}
BYTE_TO_NAME = {b: n for _, n, b, *_ in OPDEF}
BYTE_HAS_IMM = {b: imm for _, _, b, imm, *_ in OPDEF}
LETTER_HAS_IMM = {s: imm for s, _, _, imm, *_ in OPDEF}
NAME_HAS_IMM = {n: imm for _, n, _, imm, *_ in OPDEF}
ALIASES = {n: n for _, n, *_ in OPDEF}
ALIASES.update({"PUSH":"LIT","PRINT":"PUTN","EMIT":"PUTC","LOAD":"FETCH","2DUP":"DDUP","2DROP":"DDROP",">R":"TOR","R>":"RFROM"})
SYNONYM_PAIRS = [
    ("TAA","TAG","HALT"),("AAA","AAG","LIT"),("GCT","GCC","ADD"),("GAT","GAC","DUP"),
    ("GAA","GAG","EQ"),("TTT","TTC","FETCH"),("AAT","AAC","NOT"),("TGT","TGC","CALL"),
    ("TAT","TAC","PUTN"),("TCT","TCC","SUB"),("GTT","GTC","DIV"),("CGT","CGC","RET"),
]
CODON_TO_OP = {}
for a,b,name in SYNONYM_PAIRS:
    CODON_TO_OP[a] = name
    CODON_TO_OP[b] = name
CODON_TO_OP.update({
    "TTA":"JMP","TTG":"GT","TCA":"IN","TCG":"JZ","TGA":"QUIT","TGG":"PUTC",
    "CTT":"OR","CTC":"AND","CTA":"STORE","CTG":"OVER","CCT":"XOR","CCC":"SWAP",
    "CCA":"NOP","CCG":"MOD","CAT":"NEG","CAC":"DROP","CAA":"GETC","CAG":"SHL",
    "CGA":"SHR","CGG":"RAND","ATT":"INC","ATC":"DEC","ATA":"JNZ","ATG":"MUL",
    "ACT":"ROT","ACC":"NIP","ACA":"NE","ACG":"GE","AGT":"LE","AGC":"MAX",
    "AGA":"MIN","AGG":"ABS","GTA":"DEPTH","GTG":"TUCK","GCA":"PICK","GCG":"DDUP",
    "GGA":"DDROP","GGG":"TOR","GGT":"RFROM","GGC":"LT",
})
MEM_SIZE, STACK_LIMIT, CALL_LIMIT = 4096, 65536, 4096
STEP_LIMIT_DEFAULT = 1_000_000
BITS_OF_BASE = {"A":0,"C":1,"G":2,"T":3}
BASE_OF_BITS = "ACGT"

class GATCError(Exception):
    pass

@dataclass
class Token:
    kind: str
    value: object
    line: int
    col: int

@dataclass
class Instr:
    name: str
    imm: Optional[int] = None
    imm_label: Optional[str] = None
    line: int = 0

def tokenize_aa(src: str) -> List[Token]:
    tokens, i, n, line, col = [], 0, len(src), 1, 1
    def peek(k=1):
        j = i+k
        return src[j] if j < n else ""
    while i < n:
        ch = src[i]
        if ch in " \t\r":
            i += 1; col += 1; continue
        if ch == "\n":
            i += 1; line += 1; col = 1; continue
        if ch == "#":
            while i < n and src[i] != "\n":
                i += 1; col += 1
            continue
        if ch == '"':
            start = col; i += 1; col += 1; buf = []
            while i < n and src[i] != '"':
                if src[i] == "\\" and i+1 < n:
                    buf.append({"n":"\n","t":"\t","r":"\r","\\":"\\",'"':'"'}.get(src[i+1], src[i+1]))
                    i += 2; col += 2; continue
                if src[i] == "\n":
                    raise GATCError(f"line {line}:{start}: unterminated string")
                buf.append(src[i]); i += 1; col += 1
            if i >= n or src[i] != '"':
                raise GATCError(f"line {line}:{start}: unterminated string")
            i += 1; col += 1
            tokens.append(Token("STR", "".join(buf), line, start)); continue
        if ch == ":":
            start = col; i += 1; col += 1; name = []
            while i < n and (src[i].isalnum() or src[i] == "_"):
                name.append(src[i]); i += 1; col += 1
            if not name:
                raise GATCError(f"line {line}:{start}: empty label")
            tokens.append(Token("LABEL", "".join(name).upper(), line, start)); continue
        if ch.isdigit() or (ch == "-" and peek(1).isdigit()):
            start, sign = col, 1
            if ch == "-":
                sign = -1; i += 1; col += 1
            if i < n and src[i] == "0" and peek(1) in "xX":
                i += 2; col += 2; digits = []
                while i < n and src[i] in "0123456789abcdefABCDEF":
                    digits.append(src[i]); i += 1; col += 1
                num = sign * int("".join(digits), 16)
            else:
                digits = []
                while i < n and src[i].isdigit():
                    digits.append(src[i]); i += 1; col += 1
                num = sign * int("".join(digits), 10)
            tokens.append(Token("NUM", num, line, start)); continue
        if ch.isalpha() or ch in LETTER_TO_OP:
            start = col
            if ch in LETTER_TO_OP and not (peek(1).isalnum() or peek(1) == "_"):
                tokens.append(Token("OP", ch.upper() if ch.isalpha() else ch, line, start))
                i += 1; col += 1; continue
            word = []
            while i < n and (src[i].isalnum() or src[i] == "_"):
                word.append(src[i]); i += 1; col += 1
            w = "".join(word); wu = w.upper()
            if len(w) == 1 and (w.upper() if w.isalpha() else w) in LETTER_TO_OP:
                tokens.append(Token("OP", w.upper() if w.isalpha() else w, line, start))
            elif wu in ALIASES:
                tokens.append(Token("NAME", ALIASES[wu], line, start))
            else:
                tokens.append(Token("IDENT", wu, line, start))
            continue
        raise GATCError(f"line {line}:{col}: unexpected character {ch!r}")
    return tokens

def assemble_tokens(tokens):
    program, labels, i = [], {}, 0
    while i < len(tokens):
        tok = tokens[i]
        if tok.kind == "LABEL":
            if tok.value in labels:
                raise GATCError(f"line {tok.line}: duplicate label :{tok.value}")
            labels[str(tok.value)] = len(program); i += 1; continue
        if tok.kind == "STR":
            for ch in str(tok.value):
                program.append(Instr("LIT", imm=ord(ch), line=tok.line))
                program.append(Instr("PUTC", line=tok.line))
            i += 1; continue
        if tok.kind == "NUM":
            program.append(Instr("LIT", imm=int(tok.value), line=tok.line)); i += 1; continue
        if tok.kind == "IDENT":
            program.append(Instr("LIT", imm_label=str(tok.value), line=tok.line)); i += 1; continue
        if tok.kind in ("NAME","OP"):
            if tok.kind == "OP":
                name = LETTER_TO_OP[str(tok.value)]; has = LETTER_HAS_IMM[str(tok.value)]
            else:
                name = str(tok.value); has = NAME_HAS_IMM[name]
            imm = imm_label = None
            if has:
                if i+1 >= len(tokens):
                    raise GATCError(f"line {tok.line}: {name} requires an operand")
                nxt = tokens[i+1]
                if nxt.kind == "NUM": imm = int(nxt.value)
                elif nxt.kind == "IDENT": imm_label = str(nxt.value)
                else: raise GATCError(f"line {tok.line}: {name} requires number or label")
                i += 2
            else:
                i += 1
            program.append(Instr(name, imm=imm, imm_label=imm_label, line=tok.line)); continue
        raise GATCError(f"line {tok.line}: cannot assemble {tok}")
    return program, labels

def resolve(program, labels):
    out = []
    for ins in program:
        if ins.imm_label is not None:
            if ins.imm_label not in labels:
                raise GATCError(f"line {ins.line}: unknown label :{ins.imm_label}")
            out.append(Instr(ins.name, imm=labels[ins.imm_label], line=ins.line))
        else:
            out.append(ins)
    return out

def encode_bytecode(program):
    buf = bytearray()
    for ins in program:
        b = NAME_TO_BYTE[ins.name]; buf.append(b)
        if BYTE_HAS_IMM[b]:
            buf.extend((int(ins.imm) & 0xFFFFFFFF).to_bytes(4, "little"))
    return bytes(buf)

def decode_bytecode(raw):
    if raw[:4] == b"GATC":
        raw = raw[4:]
    program, i = [], 0
    while i < len(raw):
        b = raw[i]; i += 1
        name = BYTE_TO_NAME.get(b)
        if name is None:
            print(f"error: unknown opcode byte 0x{b:02X}; treating as NOP and running anyway", file=sys.stderr)
            program.append(Instr("NOP")); continue
        imm = None
        if BYTE_HAS_IMM[b]:
            if i+4 > len(raw):
                print("error: LIT truncated; pushing 0 and running anyway", file=sys.stderr)
                imm = 0; i = len(raw)
            else:
                imm = int.from_bytes(raw[i:i+4], "little"); i += 4
        program.append(Instr(name, imm=imm))
    return program

def strip_to_bases(text):
    kept = [ln for ln in text.splitlines() if not ln.startswith(">")]
    blob = "\n".join(kept) if kept else text
    return "".join(ch for ch in blob.upper().replace("U","T") if ch in "ACGT")

def bytecode_to_dna(bc):
    packet = len(bc).to_bytes(2, "big") + bc
    bases = []
    for byte in packet:
        bases.append(BASE_OF_BITS[(byte>>6)&3]); bases.append(BASE_OF_BITS[(byte>>4)&3])
        bases.append(BASE_OF_BITS[(byte>>2)&3]); bases.append(BASE_OF_BITS[byte&3])
    return "ATG" + "".join(bases) + "TAA"

def dna_to_bytecode(text):
    bases = strip_to_bases(text)
    if bases.startswith("ATG"): bases = bases[3:]
    if bases.endswith("TAA"): bases = bases[:-3]
    if len(bases) % 4:
        pad = 4 - (len(bases)%4)
        print(f"error: DNA payload is not a multiple of 4 bases; padding {pad} A and running anyway", file=sys.stderr)
        bases += "A"*pad
    raw = bytearray()
    for i in range(0, len(bases), 4):
        byte = 0
        for ch in bases[i:i+4]:
            byte = (byte<<2) | BITS_OF_BASE[ch]
        raw.append(byte)
    if len(raw) < 2:
        print("error: DNA genome missing length; running empty payload", file=sys.stderr)
        return b""
    ln = int.from_bytes(raw[:2], "big")
    payload = bytes(raw[2:])
    if ln != len(payload):
        print(f"error: DNA length header {ln} but payload is {len(payload)}; ignoring length and running anyway", file=sys.stderr)
    return payload

def codon_value(codon):
    return (BITS_OF_BASE[codon[0]]<<4)|(BITS_OF_BASE[codon[1]]<<2)|BITS_OF_BASE[codon[2]]

def compile_gene(text, require_start=True):
    bases = strip_to_bases(text)
    start = 0
    if require_start:
        idx = bases.find("ATG")
        if idx >= 0: start = idx
    seq = bases[start:]
    program, i, first = [], 0, True
    while i+2 < len(seq):
        codon = seq[i:i+3]; i += 3
        if first and codon == "ATG" and require_start:
            first = False; continue
        first = False
        name = CODON_TO_OP.get(codon)
        if name is None:
            print(f"error: unknown codon {codon}; NOP and running anyway", file=sys.stderr)
            program.append(Instr("NOP")); continue
        if name == "LIT":
            if i+2 >= len(seq):
                print("error: LIT codon missing immediate; pushing 0", file=sys.stderr)
                program.append(Instr("LIT", imm=0)); break
            imm = seq[i:i+3]; i += 3
            program.append(Instr("LIT", imm=codon_value(imm)))
        else:
            program.append(Instr(name))
        if name == "HALT":
            break
    if not program or program[-1].name != "HALT":
        program.append(Instr("HALT"))
    return program

def looks_packed(text):
    bases = strip_to_bases(text)
    if not (bases.startswith("ATG") and bases.endswith("TAA") and len(bases) >= 11):
        return False
    body = bases[3:-3]
    if len(body) < 8 or len(body) % 4:
        return False
    raw = bytearray()
    for i in range(0, 8, 4):
        byte = 0
        for ch in body[i:i+4]:
            byte = (byte<<2) | BITS_OF_BASE[ch]
        raw.append(byte)
    ln = int.from_bytes(raw, "big")
    return ln == (len(body)-8)//4

def compile_source(text, filename=""):
    lower = filename.lower()
    looks_dna = lower.endswith((".dna",".gatc",".fa",".fasta",".fna"))
    if not looks_dna:
        bases = strip_to_bases(text)
        letters = [c for c in text.upper() if c.isalpha()]
        if bases.startswith("ATG") and len(bases) >= 12 and not [c for c in letters if c not in "ACGTU"]:
            looks_dna = True
    if looks_dna:
        if looks_packed(text):
            return decode_bytecode(dna_to_bytecode(text))
        print("error: not a packed genome; running as codons anyway", file=sys.stderr)
        return compile_gene(text, require_start=False)
    program, labels = assemble_tokens(tokenize_aa(text))
    return resolve(program, labels)

def write_fasta(dna, header, width=60):
    lines = [f">{header}"]
    for i in range(0, len(dna), width):
        lines.append(dna[i:i+width])
    return "\n".join(lines)+"\n"

def translate_coding_sequence(text):
    bases = strip_to_bases(text)
    idx = bases.find("ATG")
    seq = bases[idx:] if idx >= 0 else bases
    letters, i, first = [], 0, True
    while i+2 < len(seq):
        codon = seq[i:i+3]; i += 3
        if first and codon == "ATG":
            letters.append("M"); first = False; continue
        first = False
        name = CODON_TO_OP.get(codon, "?")
        sym = next((s for s,n,*_ in OPDEF if n == name), "?")
        if name == "LIT" and i+2 < len(seq):
            imm = seq[i:i+3]; i += 3
            letters.append(f"{sym}{codon_value(imm)}")
        else:
            letters.append(sym)
        if name == "HALT":
            break
    return " ".join(letters)

def disasm(program):
    lines = []
    for pc, ins in enumerate(program):
        lines.append(f"{pc:04d}: {ins.name:<6} {ins.imm}" if ins.imm is not None else f"{pc:04d}: {ins.name}")
    return "\n".join(lines)

class VM:
    def __init__(self, program, stdin=None, stdout=None, mem_size=MEM_SIZE):
        self.program, self.pc, self.stack, self.calls = program, 0, [], []
        self.mem = [0]*mem_size
        self.stdin = stdin or sys.stdin
        self.stdout = stdout or sys.stdout
        self.halted = False
        self.exit_code = 0
        self.steps = 0
    def push(self, n):
        if len(self.stack) >= STACK_LIMIT:
            print("error: stack overflow; dropping push and running anyway", file=sys.stderr)
            return
        self.stack.append(int(n))
    def pop(self):
        if not self.stack:
            print(f"error: stack underflow at pc={self.pc}; using 0 and running anyway", file=sys.stderr)
            return 0
        return self.stack.pop()
    def run(self, max_steps=STEP_LIMIT_DEFAULT):
        n = len(self.program)
        if n == 0: return 0
        while not self.halted:
            if self.steps >= max_steps:
                print(f"error: step limit {max_steps} exceeded; stopping", file=sys.stderr)
                return self.exit_code
            if self.pc < 0 or self.pc >= n:
                print(f"error: pc {self.pc} out of range; stopping", file=sys.stderr)
                return self.exit_code
            ins = self.program[self.pc]
            self.steps += 1
            self._exec(ins)
        return self.exit_code
    def _exec(self, ins):
        op = ins.name
        if op == "LIT":
            self.push(ins.imm if ins.imm is not None else 0); self.pc += 1
        elif op == "ADD":
            b,a = self.pop(), self.pop(); self.push(a+b); self.pc += 1
        elif op == "SUB":
            b,a = self.pop(), self.pop(); self.push(a-b); self.pc += 1
        elif op == "MUL":
            b,a = self.pop(), self.pop(); self.push(a*b); self.pc += 1
        elif op == "DIV":
            b,a = self.pop(), self.pop()
            if b == 0:
                print("error: division by zero; pushing 0 and running anyway", file=sys.stderr); self.push(0)
            else:
                self.push(a//b)
            self.pc += 1
        elif op == "MOD":
            b,a = self.pop(), self.pop()
            if b == 0:
                print("error: modulo by zero; pushing 0 and running anyway", file=sys.stderr); self.push(0)
            else:
                self.push(a%b)
            self.pc += 1
        elif op == "AND":
            b,a = self.pop(), self.pop(); self.push(a&b); self.pc += 1
        elif op == "OR":
            b,a = self.pop(), self.pop(); self.push(a|b); self.pc += 1
        elif op == "XOR":
            b,a = self.pop(), self.pop(); self.push(a^b); self.pc += 1
        elif op == "NOT":
            self.push(1 if self.pop()==0 else 0); self.pc += 1
        elif op == "NEG":
            self.push(-self.pop()); self.pc += 1
        elif op == "EQ":
            b,a = self.pop(), self.pop(); self.push(1 if a==b else 0); self.pc += 1
        elif op == "LT":
            b,a = self.pop(), self.pop(); self.push(1 if a<b else 0); self.pc += 1
        elif op == "GT":
            b,a = self.pop(), self.pop(); self.push(1 if a>b else 0); self.pc += 1
        elif op == "NE":
            b,a = self.pop(), self.pop(); self.push(1 if a!=b else 0); self.pc += 1
        elif op == "GE":
            b,a = self.pop(), self.pop(); self.push(1 if a>=b else 0); self.pc += 1
        elif op == "LE":
            b,a = self.pop(), self.pop(); self.push(1 if a<=b else 0); self.pc += 1
        elif op == "MAX":
            b,a = self.pop(), self.pop(); self.push(a if a>b else b); self.pc += 1
        elif op == "MIN":
            b,a = self.pop(), self.pop(); self.push(a if a<b else b); self.pc += 1
        elif op == "ABS":
            self.push(abs(self.pop())); self.pc += 1
        elif op == "DUP":
            if not self.stack:
                print("error: stack underflow on DUP; pushing 0", file=sys.stderr); self.push(0)
            else:
                self.push(self.stack[-1])
            self.pc += 1
        elif op == "SWAP":
            b,a = self.pop(), self.pop(); self.push(b); self.push(a); self.pc += 1
        elif op == "OVER":
            if len(self.stack)<2:
                print("error: stack underflow on OVER; pushing 0", file=sys.stderr); self.push(0)
            else:
                self.push(self.stack[-2])
            self.pc += 1
        elif op == "DROP":
            self.pop(); self.pc += 1
        elif op == "NIP":
            b = self.pop(); self.pop(); self.push(b); self.pc += 1
        elif op == "ROT":
            if len(self.stack)<3:
                print("error: stack underflow on ROT; nop", file=sys.stderr)
            else:
                c,b,a = self.stack.pop(), self.stack.pop(), self.stack.pop()
                self.push(b); self.push(c); self.push(a)
            self.pc += 1
        elif op == "TUCK":
            if len(self.stack)<2:
                print("error: stack underflow on TUCK; nop", file=sys.stderr)
            else:
                self.stack.insert(-2, self.stack[-1])
            self.pc += 1
        elif op == "PICK":
            n = self.pop()
            if n < 0 or n >= len(self.stack):
                print(f"error: PICK out of range {n}; pushing 0", file=sys.stderr); self.push(0)
            else:
                self.push(self.stack[-1-n])
            self.pc += 1
        elif op == "DDUP":
            if len(self.stack)<2:
                print("error: stack underflow on DDUP; nop", file=sys.stderr)
            else:
                self.push(self.stack[-2]); self.push(self.stack[-1])
            self.pc += 1
        elif op == "DDROP":
            self.pop(); self.pop(); self.pc += 1
        elif op == "TOR":
            if len(self.calls) >= CALL_LIMIT:
                print("error: return stack overflow; dropping", file=sys.stderr); self.pop()
            else:
                self.calls.append(self.pop())
            self.pc += 1
        elif op == "RFROM":
            if not self.calls:
                print("error: return stack underflow; pushing 0 and running anyway", file=sys.stderr); self.push(0)
            else:
                self.push(self.calls.pop())
            self.pc += 1
        elif op == "INC":
            self.push(self.pop()+1); self.pc += 1
        elif op == "DEC":
            self.push(self.pop()-1); self.pc += 1
        elif op == "DEPTH":
            self.push(len(self.stack)); self.pc += 1
        elif op == "SHL":
            b,a = self.pop(), self.pop(); self.push(a << (b & 31)); self.pc += 1
        elif op == "SHR":
            b,a = self.pop(), self.pop(); self.push(a >> (b & 31)); self.pc += 1
        elif op == "FETCH":
            addr = self.pop()
            if addr < 0 or addr >= len(self.mem):
                print(f"error: fetch out of bounds {addr}; pushing 0", file=sys.stderr); self.push(0)
            else:
                self.push(self.mem[addr])
            self.pc += 1
        elif op == "STORE":
            addr, val = self.pop(), self.pop()
            if addr < 0 or addr >= len(self.mem):
                print(f"error: store out of bounds {addr}; ignoring", file=sys.stderr)
            else:
                self.mem[addr] = val
            self.pc += 1
        elif op == "PUTC":
            self.stdout.write(chr(self.pop() & 0xFF)); self.stdout.flush(); self.pc += 1
        elif op == "PUTN":
            self.stdout.write(str(self.pop())); self.stdout.flush(); self.pc += 1
        elif op == "IN":
            line = self.stdin.readline()
            if line == "": self.push(-1)
            else:
                try: self.push(int(line.strip()))
                except ValueError: self.push(0)
            self.pc += 1
        elif op == "GETC":
            ch = self.stdin.read(1); self.push(ord(ch) if ch else -1); self.pc += 1
        elif op == "RAND":
            n = self.pop(); self.push(0 if n<=0 else random.randrange(n)); self.pc += 1
        elif op == "JMP":
            self.pc = self.pop()
        elif op == "JZ":
            addr, flag = self.pop(), self.pop(); self.pc = addr if flag==0 else self.pc+1
        elif op == "JNZ":
            addr, flag = self.pop(), self.pop(); self.pc = addr if flag!=0 else self.pc+1
        elif op == "CALL":
            if len(self.calls) >= CALL_LIMIT:
                print("error: call stack overflow; ignoring call", file=sys.stderr); self.pop(); self.pc += 1
            else:
                addr = self.pop(); self.calls.append(self.pc+1); self.pc = addr
        elif op == "RET":
            if not self.calls:
                print("error: return without CALL; ignoring and running anyway", file=sys.stderr)
                self.pc += 1
            else:
                self.pc = self.calls.pop()
        elif op == "NOP":
            self.pc += 1
        elif op == "HALT":
            self.halted = True; self.exit_code = 0
        elif op == "QUIT":
            self.exit_code = self.pop(); self.halted = True
        else:
            print(f"error: unimplemented opcode {op}; nop", file=sys.stderr); self.pc += 1

def cmd_opcodes():
    by = {}
    for codon, name in sorted(CODON_TO_OP.items()):
        by.setdefault(name, []).append(codon)
    print(f"HolyGATC opcode table — {len(OPDEF)} opcodes, {len(CODON_TO_OP)} codons\n")
    print(f"{'Sym':<4} {'Name':<7} {'Byte':<6} {'Imm':<5} {'Codons':<12} Stack")
    print("-"*78)
    for letter, name, byte, imm, stack, _ in OPDEF:
        print(f"{letter:<4} {name:<7} 0x{byte:02X}  {str(imm):<5} {','.join(by.get(name,['-'])):<12} {stack}")

def main(argv=None):
    p = argparse.ArgumentParser(prog="gatc52", description="HolyGATC 52 — 52 opcodes, 64 codons.")
    p.add_argument("command", nargs="?", default="help", choices=["run","compile","disasm","translate","opcodes","help","version"])
    p.add_argument("source", nargs="?", help="source file (.aa / .dna / .fa / .fasta)")
    p.add_argument("-o","--output")
    p.add_argument("--dna", action="store_true")
    p.add_argument("--max-steps", type=int, default=STEP_LIMIT_DEFAULT)
    p.add_argument("--no-banner", action="store_true")
    p.add_argument("--gene", action="store_true")
    args = p.parse_intermixed_args(argv)
    if args.command == "help" and not args.source:
        print(BANNER); p.print_help()
        print("\nExamples:\n  python3 GATC52.py run program.aa --no-banner\n  python3 GATC52.py compile program.aa -o program.dna --dna\n  python3 GATC52.py run human.dna --no-banner\n  python3 GATC52.py opcodes")
        return 0
    if args.command == "version":
        print(f"HolyGATC {VERSION} — 52 opcodes — 64 codons"); return 0
    if args.command == "opcodes":
        cmd_opcodes(); return 0
    if not args.source:
        print("error: source file required", file=sys.stderr); return 2
    try:
        with open(args.source, encoding="utf-8", errors="replace") as f:
            text = f.read()
        if args.command == "translate":
            sys.stdout.write(translate_coding_sequence(text)+"\n"); return 0
        if args.command == "compile":
            program = compile_gene(text) if args.gene else compile_source(text, args.source)
            bc = encode_bytecode(program)
            if args.dna:
                payload = write_fasta(bytecode_to_dna(bc), f"HolyGATC52 genome from {args.source}")
                if args.output:
                    open(args.output,"w",encoding="utf-8").write(payload)
                else:
                    sys.stdout.write(payload)
            else:
                if args.output:
                    open(args.output,"wb").write(bc)
                else:
                    sys.stdout.buffer.write(bc)
            return 0
        if args.command == "disasm":
            if args.gene: program = compile_gene(text)
            elif args.source.lower().endswith(".gbc"): program = decode_bytecode(open(args.source,"rb").read())
            else: program = compile_source(text, args.source)
            sys.stdout.write(disasm(program)+"\n"); return 0
        if not args.no_banner:
            print(BANNER, file=sys.stderr)
        if args.gene: program = compile_gene(text, require_start=False)
        elif args.source.lower().endswith(".gbc"): program = decode_bytecode(open(args.source,"rb").read())
        else: program = compile_source(text, args.source)
        return VM(program).run(max_steps=args.max_steps)
    except GATCError as e:
        print(f"error: {e}", file=sys.stderr); return 1

if __name__ == "__main__":
    sys.exit(main())
