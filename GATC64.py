#!/usr/bin/env python3
"""HolyGATC / GATC64 version 4.0.0
64 opcodes, one per codon. Real FASTA runs anyway.
Packed DNA: ATG + 16-bit big-endian length + payload + TAA.
"""
from __future__ import annotations
import argparse, random, sys
from dataclasses import dataclass
from typing import Optional
VERSION = "4.0.0"
BANNER = r"""
  ____    _  _____ ____  __ _  _
 / ___|  / \|_   _/ ___|/ /| || |
| |  _  / _ \ | || |   / / | || |_
| |_| |/ ___ \| || |__/ /  |__   _|
 \____/_/   \_\_| \____/      |_|
  64 opcodes. 64 codons. One each.
  Version 4. Real FASTA runs anyway.
"""
OPDEF = [
 ("A","ADD",0x01,False,"(a b -- a+b)"),("B","JMP",0x02,False,"(addr --)"),
 ("C","CALL",0x03,False,"(addr --)"),("D","DUP",0x04,False,"(a -- a a)"),
 ("E","EQ",0x05,False,"(a b -- flag)"),("F","FETCH",0x06,False,"(addr -- val)"),
 ("G","GT",0x07,False,"(a b -- flag)"),("H","HALT",0x08,False,"(--)"),
 ("I","IN",0x09,False,"(-- n)"),("J","JZ",0x0A,False,"(flag addr --)"),
 ("K","LIT",0x0B,True,"(-- n)"),("L","LT",0x0C,False,"(a b -- flag)"),
 ("M","MUL",0x0D,False,"(a b -- a*b)"),("N","NOT",0x0E,False,"(a -- flag)"),
 ("O","OR",0x0F,False,"(a b -- a|b)"),("P","AND",0x10,False,"(a b -- a&b)"),
 ("Q","QUIT",0x11,False,"(code --)"),("R","RET",0x12,False,"(--)"),
 ("S","SUB",0x13,False,"(a b -- a-b)"),("T","STORE",0x14,False,"(val addr --)"),
 ("U","OVER",0x15,False,"(a b -- a b a)"),("V","DIV",0x16,False,"(a b -- a/b)"),
 ("W","PUTC",0x17,False,"(n --)"),("X","XOR",0x18,False,"(a b -- a^b)"),
 ("Y","PUTN",0x19,False,"(n --)"),("Z","SWAP",0x1A,False,"(a b -- b a)"),
 ("_","NOP",0x1B,False,"(--)"),("%","MOD",0x1C,False,"(a b -- a%b)"),
 ("~","NEG",0x1D,False,"(a -- -a)"),("!","DROP",0x1E,False,"(a --)"),
 ("?","GETC",0x1F,False,"(-- n)"),("<","SHL",0x20,False,"(a b -- a<<b)"),
 (">","SHR",0x21,False,"(a b -- a>>b)"),("$","RAND",0x22,False,"(n -- r)"),
 ("+","INC",0x23,False,"(a -- a+1)"),("-","DEC",0x24,False,"(a -- a-1)"),
 ("*","JNZ",0x25,False,"(flag addr --)"),("@","ROT",0x26,False,"(a b c -- b c a)"),
 ("&","NIP",0x27,False,"(a b -- b)"),("=","NE",0x28,False,"(a b -- flag)"),
 ("/","GE",0x29,False,"(a b -- flag)"),("^","LE",0x2A,False,"(a b -- flag)"),
 ("|","MAX",0x2B,False,"(a b -- max)"),(";","MIN",0x2C,False,"(a b -- min)"),
 (".","ABS",0x2D,False,"(a -- |a|)"),(",","DEPTH",0x2E,False,"(-- n)"),
 ("`","TUCK",0x2F,False,"(a b -- b a b)"),("'","PICK",0x30,False,"(n -- x)"),
 ("[","DDUP",0x31,False,"(a b -- a b a b)"),("]","DDROP",0x32,False,"(a b --)"),
 ("{","TOR",0x33,False,"(a --)"),("}","RFROM",0x34,False,"(-- a)"),
 ("\\","INV",0x35,False,"(a -- ~a)"),("(","ZERO",0x36,False,"(-- 0)"),
 (")","ONE",0x37,False,"(-- 1)"),("CR","CR",0x38,False,"(--) newline"),
 ("TSWP","TSWP",0x39,False,"(a b c d -- c d a b)"),("BOOL","BOOL",0x3A,False,"(a -- 0|1)"),
 ("PSTO","PSTO",0x3B,False,"(addr --) mem++"),("UNDER","UNDER",0x3C,False,"(a b -- a a b)"),
 ("RDROP","RDROP",0x3D,False,"(R: a --)"),("NROLL","NROLL",0x3E,False,"(a b c -- c a b)"),
 ("SIGN","SIGN",0x3F,False,"(a -- -1|0|1)"),("MADD","MADD",0x40,False,"(n addr --)"),
]
LETTER_TO_OP={s:n for s,n,*_ in OPDEF if len(s)==1}
NAME_TO_BYTE={n:b for _,n,b,*_ in OPDEF}
BYTE_TO_NAME={b:n for _,n,b,*_ in OPDEF}
BYTE_HAS_IMM={b:i for _,_,b,i,*_ in OPDEF}
LETTER_HAS_IMM={s:i for s,_,_,i,*_ in OPDEF if len(s)==1}
NAME_HAS_IMM={n:i for _,n,_,i,*_ in OPDEF}
ALIASES={n:n for _,n,*_ in OPDEF}
ALIASES.update({"PUSH":"LIT","PRINT":"PUTN","EMIT":"PUTC","LOAD":"FETCH","2DUP":"DDUP","2DROP":"DDROP",">R":"TOR","R>":"RFROM","2SWAP":"TSWP","INVERT":"INV"})
CODON_TO_OP={
"TTT":"FETCH","TTC":"LT","TTA":"JMP","TTG":"GT","TCT":"SUB","TCC":"IN","TCA":"JZ","TCG":"PUTC",
"TAT":"PUTN","TAC":"QUIT","TAA":"HALT","TAG":"NOP","TGT":"CALL","TGC":"RET","TGA":"OR","TGG":"AND",
"CTT":"STORE","CTC":"OVER","CTA":"XOR","CTG":"SWAP","CCT":"MOD","CCC":"NEG","CCA":"DROP","CCG":"GETC",
"CAT":"SHL","CAC":"SHR","CAA":"RAND","CAG":"INC","CGT":"DEC","CGC":"JNZ","CGA":"ROT","CGG":"NIP",
"ATT":"NE","ATC":"GE","ATA":"LE","ATG":"MUL","ACT":"MAX","ACC":"MIN","ACA":"ABS","ACG":"DEPTH",
"AAT":"NOT","AAC":"EQ","AAA":"LIT","AAG":"ADD","AGT":"DUP","AGC":"DIV","AGA":"TUCK","AGG":"PICK",
"GTT":"DDUP","GTC":"DDROP","GTA":"TOR","GTG":"RFROM","GCT":"INV","GCC":"ZERO","GCA":"ONE","GCG":"CR",
"GAT":"TSWP","GAC":"BOOL","GAA":"PSTO","GAG":"UNDER","GGT":"RDROP","GGC":"NROLL","GGA":"SIGN","GGG":"MADD"}
MEM_SIZE,STACK_LIMIT,CALL_LIMIT,STEP_LIMIT_DEFAULT=4096,65536,4096,5_000_000
BITS={"A":0,"C":1,"G":2,"T":3}; BASE="ACGT"
class GATCError(Exception): pass
_seen=set()
def warn(msg):
    key=msg.split(";")[0][:48]
    if key in _seen: return
    _seen.add(key); print("error: "+msg, file=sys.stderr)
@dataclass
class Token:
    kind:str; value:object; line:int; col:int
@dataclass
class Instr:
    name:str; imm:Optional[int]=None; imm_label:Optional[str]=None; line:int=0
def tokenize_aa(src):
    tokens,i,n,line,col=[],0,len(src),1,1
    def peek(k=1):
        j=i+k; return src[j] if j<n else ""
    while i<n:
        ch=src[i]
        if ch in " \t\r":
            i+=1; col+=1; continue
        if ch=="\n":
            i+=1; line+=1; col=1; continue
        if ch=="#":
            while i<n and src[i]!="\n": i+=1; col+=1
            continue
        if ch=='"':
            start=col; i+=1; col+=1; buf=[]
            while i<n and src[i]!='"':
                if src[i]=="\\" and i+1<n:
                    buf.append({"n":"\n","t":"\t","r":"\r","\\":"\\",'"':'"'}.get(src[i+1],src[i+1])); i+=2; col+=2; continue
                if src[i]=="\n": raise GATCError(f"line {line}: unterminated string")
                buf.append(src[i]); i+=1; col+=1
            if i>=n or src[i]!='"': raise GATCError("unterminated string")
            i+=1; col+=1; tokens.append(Token("STR","".join(buf),line,start)); continue
        if ch==":":
            start=col; i+=1; col+=1; name=[]
            while i<n and (src[i].isalnum() or src[i]=="_"): name.append(src[i]); i+=1; col+=1
            if not name: raise GATCError("empty label")
            tokens.append(Token("LABEL","".join(name).upper(),line,start)); continue
        if ch.isdigit() or (ch=="-" and peek(1).isdigit()):
            start,sign=col,1
            if ch=="-": sign=-1; i+=1; col+=1
            if i<n and src[i]=="0" and peek(1) in "xX":
                i+=2; col+=2; d=[]
                while i<n and src[i] in "0123456789abcdefABCDEF": d.append(src[i]); i+=1; col+=1
                num=sign*int("".join(d),16)
            else:
                d=[]
                while i<n and src[i].isdigit(): d.append(src[i]); i+=1; col+=1
                num=sign*int("".join(d) or "0",10)
            tokens.append(Token("NUM",num,line,start)); continue
        if ch.isalpha() or ch in LETTER_TO_OP:
            start=col
            if ch in LETTER_TO_OP and not (peek(1).isalnum() or peek(1)=="_"):
                tokens.append(Token("OP", ch.upper() if ch.isalpha() else ch, line, start)); i+=1; col+=1; continue
            w=[]
            while i<n and (src[i].isalnum() or src[i]=="_"): w.append(src[i]); i+=1; col+=1
            word="".join(w); wu=word.upper()
            if len(word)==1 and (word.upper() if word.isalpha() else word) in LETTER_TO_OP:
                tokens.append(Token("OP", word.upper() if word.isalpha() else word, line, start))
            elif wu in ALIASES: tokens.append(Token("NAME", ALIASES[wu], line, start))
            else: tokens.append(Token("IDENT", wu, line, start))
            continue
        raise GATCError(f"line {line}:{col}: unexpected {ch!r}")
    return tokens
def assemble_tokens(tokens):
    program,labels,i=[],{},0
    while i<len(tokens):
        tok=tokens[i]
        if tok.kind=="LABEL":
            if tok.value in labels: raise GATCError(f"duplicate :{tok.value}")
            labels[str(tok.value)]=len(program); i+=1; continue
        if tok.kind=="STR":
            for ch in str(tok.value):
                program += [Instr("LIT",imm=ord(ch),line=tok.line), Instr("PUTC",line=tok.line)]
            i+=1; continue
        if tok.kind=="NUM":
            program.append(Instr("LIT",imm=int(tok.value),line=tok.line)); i+=1; continue
        if tok.kind=="IDENT":
            program.append(Instr("LIT",imm_label=str(tok.value),line=tok.line)); i+=1; continue
        if tok.kind in ("NAME","OP"):
            name = LETTER_TO_OP[str(tok.value)] if tok.kind=="OP" else str(tok.value)
            has = LETTER_HAS_IMM[str(tok.value)] if tok.kind=="OP" else NAME_HAS_IMM[name]
            imm=imm_label=None
            if has:
                if i+1>=len(tokens): raise GATCError(f"{name} needs an operand")
                nxt=tokens[i+1]
                if nxt.kind=="NUM": imm=int(nxt.value)
                elif nxt.kind=="IDENT": imm_label=str(nxt.value)
                else: raise GATCError(f"{name} needs number or label")
                i+=2
            else: i+=1
            program.append(Instr(name,imm=imm,imm_label=imm_label,line=tok.line)); continue
        raise GATCError(f"cannot assemble {tok}")
    return program, labels
def resolve(program, labels):
    out=[]
    for ins in program:
        if ins.imm_label is not None:
            if ins.imm_label not in labels: raise GATCError(f"unknown label :{ins.imm_label}")
            out.append(Instr(ins.name, imm=labels[ins.imm_label], line=ins.line))
        else: out.append(ins)
    return out
def encode_bytecode(program):
    buf=bytearray()
    for ins in program:
        b=NAME_TO_BYTE[ins.name]; buf.append(b)
        if BYTE_HAS_IMM[b]: buf.extend((int(ins.imm or 0)&0xFFFFFFFF).to_bytes(4,"little"))
    return bytes(buf)
def decode_bytecode(raw):
    if raw[:4]==b"GATC": raw=raw[4:]
    program,i=[],0
    while i<len(raw):
        b=raw[i]; i+=1
        name=BYTE_TO_NAME.get(b)
        if name is None:
            warn(f"unknown opcode byte 0x{b:02X}; NOP"); program.append(Instr("NOP")); continue
        imm=None
        if BYTE_HAS_IMM[b]:
            if i+4>len(raw): warn("LIT truncated; pushing 0"); imm=0; i=len(raw)
            else: imm=int.from_bytes(raw[i:i+4],"little"); i+=4
        program.append(Instr(name, imm=imm))
    return program
def strip_to_bases(text):
    kept=[ln for ln in text.splitlines() if not ln.startswith(">")]
    blob="\n".join(kept) if kept else text
    return "".join(ch for ch in blob.upper().replace("U","T") if ch in "ACGT")
def bytecode_to_dna(bc):
    packet=len(bc).to_bytes(2,"big")+bc
    bases=[]
    for byte in packet:
        bases += [BASE[(byte>>6)&3],BASE[(byte>>4)&3],BASE[(byte>>2)&3],BASE[byte&3]]
    return "ATG"+"".join(bases)+"TAA"
def dna_to_bytecode(text):
    bases=strip_to_bases(text)
    if bases.startswith("ATG"): bases=bases[3:]
    if bases.endswith("TAA"): bases=bases[:-3]
    if len(bases)%4:
        pad=4-len(bases)%4; warn(f"payload not multiple of 4; padding {pad} A"); bases+="A"*pad
    raw=bytearray()
    for i in range(0,len(bases),4):
        byte=0
        for ch in bases[i:i+4]: byte=(byte<<2)|BITS[ch]
        raw.append(byte)
    if len(raw)<2: warn("missing length; empty payload"); return b""
    ln=int.from_bytes(raw[:2],"big"); payload=bytes(raw[2:])
    if ln!=len(payload): warn(f"length header {ln} but payload is {len(payload)}; ignoring length")
    return payload
def codon_value(c):
    return (BITS[c[0]]<<4)|(BITS[c[1]]<<2)|BITS[c[2]]
def compile_gene(text, force_start=True):
    bases=strip_to_bases(text)
    if force_start and not bases.startswith("ATG"):
        warn("no start codon; forcing ATG"); bases="ATG"+bases
    program,i,first=[],0,True
    while i+2<len(bases):
        codon=bases[i:i+3]; i+=3
        if first and codon=="ATG":
            first=False; continue
        first=False
        name=CODON_TO_OP.get(codon,"NOP")
        if name=="LIT":
            if i+2>=len(bases):
                warn("LIT missing immediate; pushing 0"); program.append(Instr("LIT",imm=0)); break
            imm=bases[i:i+3]; i+=3; program.append(Instr("LIT", imm=codon_value(imm)))
        else: program.append(Instr(name))
    if not program: program.append(Instr("HALT"))
    return program
def looks_packed(text):
    bases=strip_to_bases(text)
    if not (bases.startswith("ATG") and bases.endswith("TAA") and len(bases)>=11): return False
    body=bases[3:-3]
    if len(body)<8 or len(body)%4: return False
    raw=bytearray()
    for i in range(0,8,4):
        byte=0
        for ch in body[i:i+4]: byte=(byte<<2)|BITS[ch]
        raw.append(byte)
    return int.from_bytes(raw,"big")== (len(body)-8)//4
def compile_source(text, filename=""):
    lower=filename.lower()
    dna=lower.endswith((".dna",".gatc",".fa",".fasta",".fna"))
    if not dna:
        bases=strip_to_bases(text)
        letters=[c for c in text.upper() if c.isalpha()]
        if bases.startswith("ATG") and len(bases)>=12 and not [c for c in letters if c not in "ACGTU"]:
            dna=True
    if dna:
        if looks_packed(text): return decode_bytecode(dna_to_bytecode(text)), False
        warn("not a packed genome; running every codon")
        return compile_gene(text, force_start=True), True
    program, labels = assemble_tokens(tokenize_aa(text))
    return resolve(program, labels), False
def write_fasta(dna, header, width=60):
    lines=[">"+header]
    for i in range(0,len(dna),width): lines.append(dna[i:i+width])
    return "\n".join(lines)+"\n"
def translate_coding_sequence(text):
    bases=strip_to_bases(text); forced=False
    if not bases.startswith("ATG"): bases="ATG"+bases; forced=True
    out,i,first=[],0,True
    while i+2<len(bases):
        codon=bases[i:i+3]; i+=3
        if first and codon=="ATG":
            out.append("M"); first=False; continue
        first=False
        name=CODON_TO_OP.get(codon,"?")
        sym=next((s for s,n,*_ in OPDEF if n==name), name[:1])
        if name=="LIT" and i+2<len(bases):
            imm=bases[i:i+3]; i+=3; out.append(f"K{codon_value(imm)}")
        else: out.append(sym if len(sym)==1 else name)
    return ("FORCED " if forced else "")+" ".join(out)
def disasm(program):
    return "\n".join(f"{pc:04d}: {ins.name:<6} {ins.imm}" if ins.imm is not None else f"{pc:04d}: {ins.name}" for pc,ins in enumerate(program))
class VM:
    def __init__(self, program, genome=False, stdin=None, stdout=None):
        self.program=program; self.pc=0; self.stack=[]; self.calls=[]; self.mem=[0]*MEM_SIZE
        self.stdin=stdin or sys.stdin; self.stdout=stdout or sys.stdout
        self.halted=False; self.exit_code=0; self.steps=0; self.genome=genome; self._halt_note=False
    def push(self,n):
        if len(self.stack)>=STACK_LIMIT: warn("stack overflow; dropping push"); return
        self.stack.append(int(n))
    def pop(self):
        if not self.stack: warn(f"stack underflow at pc={self.pc}; using 0"); return 0
        return self.stack.pop()
    def run(self, max_steps=STEP_LIMIT_DEFAULT):
        n=len(self.program)
        if n==0: return 0
        while not self.halted:
            if self.steps>=max_steps: warn(f"step limit {max_steps} exceeded; stopping"); return self.exit_code
            if self.pc<0 or self.pc>=n: warn(f"pc {self.pc} out of range; stopping"); return self.exit_code
            self.steps+=1; self._exec(self.program[self.pc])
        return self.exit_code
    def _exec(self, ins):
        op=ins.name
        if op=="LIT": self.push(0 if ins.imm is None else ins.imm); self.pc+=1
        elif op=="ADD":
            b,a=self.pop(),self.pop(); self.push(a+b); self.pc+=1
        elif op=="SUB":
            b,a=self.pop(),self.pop(); self.push(a-b); self.pc+=1
        elif op=="MUL":
            b,a=self.pop(),self.pop(); self.push(a*b); self.pc+=1
        elif op=="DIV":
            b,a=self.pop(),self.pop(); self.push(0 if b==0 else a//b)
            if b==0: warn("division by zero; pushing 0")
            self.pc+=1
        elif op=="MOD":
            b,a=self.pop(),self.pop(); self.push(0 if b==0 else a%b)
            if b==0: warn("modulo by zero; pushing 0")
            self.pc+=1
        elif op in ("AND","OR","XOR"):
            b,a=self.pop(),self.pop(); self.push({"AND":a&b,"OR":a|b,"XOR":a^b}[op]); self.pc+=1
        elif op=="INV": self.push(~self.pop()); self.pc+=1
        elif op=="NOT": self.push(1 if self.pop()==0 else 0); self.pc+=1
        elif op=="NEG": self.push(-self.pop()); self.pc+=1
        elif op=="EQ":
            b,a=self.pop(),self.pop(); self.push(1 if a==b else 0); self.pc+=1
        elif op=="LT":
            b,a=self.pop(),self.pop(); self.push(1 if a<b else 0); self.pc+=1
        elif op=="GT":
            b,a=self.pop(),self.pop(); self.push(1 if a>b else 0); self.pc+=1
        elif op=="NE":
            b,a=self.pop(),self.pop(); self.push(1 if a!=b else 0); self.pc+=1
        elif op=="GE":
            b,a=self.pop(),self.pop(); self.push(1 if a>=b else 0); self.pc+=1
        elif op=="LE":
            b,a=self.pop(),self.pop(); self.push(1 if a<=b else 0); self.pc+=1
        elif op=="MAX":
            b,a=self.pop(),self.pop(); self.push(a if a>b else b); self.pc+=1
        elif op=="MIN":
            b,a=self.pop(),self.pop(); self.push(a if a<b else b); self.pc+=1
        elif op=="ABS": self.push(abs(self.pop())); self.pc+=1
        elif op=="SIGN":
            a=self.pop(); self.push((a>0)-(a<0)); self.pc+=1
        elif op=="BOOL": self.push(0 if self.pop()==0 else 1); self.pc+=1
        elif op=="ZERO": self.push(0); self.pc+=1
        elif op=="ONE": self.push(1); self.pc+=1
        elif op=="DUP": self.push(self.stack[-1] if self.stack else 0); self.pc+=1
        elif op=="SWAP":
            b,a=self.pop(),self.pop(); self.push(b); self.push(a); self.pc+=1
        elif op=="OVER": self.push(self.stack[-2] if len(self.stack)>=2 else 0); self.pc+=1
        elif op=="DROP": self.pop(); self.pc+=1
        elif op=="NIP":
            b=self.pop(); self.pop(); self.push(b); self.pc+=1
        elif op=="ROT":
            if len(self.stack)<3: warn("ROT underflow")
            else:
                c,b,a=self.stack.pop(),self.stack.pop(),self.stack.pop(); self.push(b); self.push(c); self.push(a)
            self.pc+=1
        elif op=="NROLL":
            if len(self.stack)<3: warn("NROLL underflow")
            else:
                c,b,a=self.stack.pop(),self.stack.pop(),self.stack.pop(); self.push(c); self.push(a); self.push(b)
            self.pc+=1
        elif op=="TUCK":
            if len(self.stack)<2: warn("TUCK underflow")
            else: self.stack.insert(-2, self.stack[-1])
            self.pc+=1
        elif op=="UNDER":
            if len(self.stack)<2: warn("UNDER underflow")
            else: self.stack.insert(-1, self.stack[-2])
            self.pc+=1
        elif op=="PICK":
            n=self.pop(); self.push(self.stack[-1-n] if 0<=n<len(self.stack) else 0); self.pc+=1
        elif op=="DDUP":
            if len(self.stack)<2: warn("DDUP underflow")
            else: self.push(self.stack[-2]); self.push(self.stack[-1])
            self.pc+=1
        elif op=="DDROP": self.pop(); self.pop(); self.pc+=1
        elif op=="TSWP":
            if len(self.stack)<4: warn("TSWP underflow")
            else:
                d,c,b,a=self.stack.pop(),self.stack.pop(),self.stack.pop(),self.stack.pop()
                self.push(c); self.push(d); self.push(a); self.push(b)
            self.pc+=1
        elif op=="TOR":
            if len(self.calls)>=CALL_LIMIT: warn("return stack overflow"); self.pop()
            else: self.calls.append(self.pop())
            self.pc+=1
        elif op=="RFROM":
            self.push(self.calls.pop() if self.calls else 0); self.pc+=1
        elif op=="RDROP":
            if self.calls: self.calls.pop()
            else: warn("RDROP on empty return stack")
            self.pc+=1
        elif op=="INC": self.push(self.pop()+1); self.pc+=1
        elif op=="DEC": self.push(self.pop()-1); self.pc+=1
        elif op=="DEPTH": self.push(len(self.stack)); self.pc+=1
        elif op=="SHL":
            b,a=self.pop(),self.pop(); self.push(a<<(b&31)); self.pc+=1
        elif op=="SHR":
            b,a=self.pop(),self.pop(); self.push(a>>(b&31)); self.pc+=1
        elif op=="FETCH":
            addr=self.pop(); self.push(self.mem[addr] if 0<=addr<len(self.mem) else 0); self.pc+=1
        elif op=="STORE":
            addr,val=self.pop(),self.pop()
            if 0<=addr<len(self.mem): self.mem[addr]=val
            else: warn(f"store out of bounds {addr}")
            self.pc+=1
        elif op=="PSTO":
            addr=self.pop()
            if 0<=addr<len(self.mem): self.mem[addr]+=1
            else: warn(f"psto out of bounds {addr}")
            self.pc+=1
        elif op=="MADD":
            addr,n=self.pop(),self.pop()
            if 0<=addr<len(self.mem): self.mem[addr]+=n
            else: warn(f"madd out of bounds {addr}")
            self.pc+=1
        elif op=="PUTC":
            self.stdout.write(chr(self.pop()&0xFF)); self.stdout.flush(); self.pc+=1
        elif op=="PUTN":
            self.stdout.write(str(self.pop())); self.stdout.flush(); self.pc+=1
        elif op=="CR":
            self.stdout.write("\n"); self.stdout.flush(); self.pc+=1
        elif op=="IN":
            line=self.stdin.readline()
            if line=="": self.push(-1)
            else:
                try: self.push(int(line.strip()))
                except ValueError: self.push(0)
            self.pc+=1
        elif op=="GETC":
            ch=self.stdin.read(1); self.push(ord(ch) if ch else -1); self.pc+=1
        elif op=="RAND":
            n=self.pop(); self.push(0 if n<=0 else random.randrange(n)); self.pc+=1
        elif op=="JMP":
            addr=self.pop()
            if self.genome: warn("JMP in genome; walking forward"); self.pc+=1
            elif isinstance(addr,int) and 0<=addr<len(self.program): self.pc=addr
            else: warn(f"bad jump {addr}; skipping"); self.pc+=1
        elif op=="JZ":
            addr,flag=self.pop(),self.pop()
            if self.genome: warn("JZ in genome; walking forward"); self.pc+=1
            else: self.pc = addr if flag==0 and 0<=addr<len(self.program) else self.pc+1
        elif op=="JNZ":
            addr,flag=self.pop(),self.pop()
            if self.genome: warn("JNZ in genome; walking forward"); self.pc+=1
            else: self.pc = addr if flag!=0 and 0<=addr<len(self.program) else self.pc+1
        elif op=="CALL":
            addr=self.pop()
            if self.genome: warn("CALL in genome; walking forward"); self.pc+=1
            elif len(self.calls)>=CALL_LIMIT or not (0<=addr<len(self.program)):
                warn(f"bad call {addr}; skipping"); self.pc+=1
            else: self.calls.append(self.pc+1); self.pc=addr
        elif op=="RET":
            if self.genome or not self.calls:
                if not self.calls: warn("return without CALL; ignoring")
                self.pc+=1
            else: self.pc=self.calls.pop()
        elif op=="NOP": self.pc+=1
        elif op=="HALT":
            if self.genome:
                if not self._halt_note: warn("HALT in genome; continuing"); self._halt_note=True
                self.pc+=1
            else: self.halted=True
        elif op=="QUIT":
            code=self.pop()
            if self.genome: warn(f"QUIT {code} in genome; continuing"); self.pc+=1
            else: self.exit_code=code; self.halted=True
        else: warn(f"unknown op {op}; nop"); self.pc+=1
def cmd_opcodes():
    by={}
    for c,n in sorted(CODON_TO_OP.items()): by.setdefault(n,[]).append(c)
    print(f"HolyGATC {VERSION} — {len(OPDEF)} opcodes, {len(CODON_TO_OP)} codons\n")
    print(f"{'Sym':<6} {'Name':<7} {'Byte':<6} {'Codon':<8} Stack")
    print("-"*62)
    for sym,name,byte,imm,stack in OPDEF:
        print(f"{sym:<6} {name:<7} 0x{byte:02X}  {','.join(by.get(name,['-'])):<8} {stack}")
def main(argv=None):
    p=argparse.ArgumentParser(prog="gatc64")
    p.add_argument("command", nargs="?", default="help", choices=["run","compile","disasm","translate","opcodes","help","version"])
    p.add_argument("source", nargs="?")
    p.add_argument("-o","--output"); p.add_argument("--dna", action="store_true")
    p.add_argument("--max-steps", type=int, default=STEP_LIMIT_DEFAULT)
    p.add_argument("--no-banner", action="store_true"); p.add_argument("--gene", action="store_true")
    args=p.parse_intermixed_args(argv)
    if args.command=="help" and not args.source:
        print(BANNER); p.print_help()
        print("\n  python3 GATC64.py run program.aa --no-banner")
        print("  python3 GATC64.py compile program.aa --dna -o program.dna")
        print("  python3 GATC64.py run human.dna --no-banner"); return 0
    if args.command=="version":
        print(f"HolyGATC {VERSION} — 64 opcodes — 64 codons"); return 0
    if args.command=="opcodes":
        cmd_opcodes(); return 0
    if not args.source: warn("source file required"); return 2
    try:
        text=open(args.source, encoding="utf-8", errors="replace").read()
        if args.command=="translate":
            sys.stdout.write(translate_coding_sequence(text)+"\n"); return 0
        if args.command=="compile":
            program,_ = (compile_gene(text), True) if args.gene else compile_source(text, args.source)
            bc=encode_bytecode(program)
            if args.dna:
                payload=write_fasta(bytecode_to_dna(bc), f"HolyGATC64 genome from {args.source}")
                (open(args.output,"w",encoding="utf-8").write(payload) if args.output else sys.stdout.write(payload))
            else:
                (open(args.output,"wb").write(bc) if args.output else sys.stdout.buffer.write(bc))
            return 0
        if args.command=="disasm":
            if args.gene: program=compile_gene(text)
            elif args.source.lower().endswith(".gbc"): program=decode_bytecode(open(args.source,"rb").read())
            else: program,_=compile_source(text, args.source)
            sys.stdout.write(disasm(program)+"\n"); return 0
        if not args.no_banner: print(BANNER, file=sys.stderr)
        if args.gene: program, genome = compile_gene(text), True
        elif args.source.lower().endswith(".gbc"): program, genome = decode_bytecode(open(args.source,"rb").read()), False
        else: program, genome = compile_source(text, args.source)
        return VM(program, genome=genome).run(max_steps=args.max_steps)
    except GATCError as e:
        warn(str(e)); return 1
if __name__=="__main__":
    sys.exit(main())
