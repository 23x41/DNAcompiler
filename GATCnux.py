#!/usr/bin/env python3
"""HolyGATC / GATC64 version 4.0.0 + file load support
64 opcodes, one per codon. Real FASTA runs anyway.
Packed DNA: ATG + 16-bit big-endian length + payload + TAA.
"""
from __future__ import annotations
import argparse, random, sys
from dataclasses import dataclass
from typing import Optional
import os, zlib, collections
import socket, ssl, struct, time, ipaddress, re
VERSION = "4.2.0-fs"
BANNER = r"""
  ____    _  _____ ____  __ _  _
 / ___|  / \|_   _/ ___|/ /| || |
| |  _  / _ \ | || |   / / | || |_
| |_| |/ ___ \| || |__/ /  |__   _|
 \____/_/   \_\_| \____/      |_|
  64 opcodes. 64 codons. One each.
  Version 4 + --load/--file/--net/--fs. Real FASTA runs anyway.
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
# --- file-backed address space (used by hexeditor.aa) -----------------------------------
# FETCH/STORE at FILE_BASE+n read/write byte n of the file given with --file (up to 1 GiB).
# Cells 4096.. are control cells.  STORE a command number into CTL_CMD to run a bulk operation
# in the host (fast); arguments go in CTL_A1/A2/A3, results come back in CTL_R1/R2/R3.
FILE_BASE=65536
MAX_FILE=1<<30                 # 1 GiB editing limit
GROW_CAP=1<<28                 # largest single insert (256 MiB)
UNDO_CAP=1<<28                 # undo history budget (about 256 MiB)
CHUNK=1<<24
CTL_CMD,CTL_A1,CTL_A2,CTL_SIZE,CTL_DIRTY,CTL_MAGIC,CTL_STATUS,CTL_R1,CTL_R2,CTL_R3,CTL_A3,CTL_LOST=range(4096,4108)
# commands: 1 save  2 insert zeros(A1=off,A2=n)  3 delete(A1,A2)  4 find(A1=start,A2=len,A3=mem addr of pattern)
#           5 checksum(A1,A2 -> R1 sum8,R2 xor8,R3 crc32)  6 undo  7 begin undo group
#           8 fill(A1,A2,A3=value)  9 copy(A1=src,A2=len,A3=dst)  10 xor(A1,A2,A3=key)  11 set size(A1)
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
# ---------------------------------------------------------------------------------------
# networking for netkit.aa  (opt-in: run with --net)
#
# The VM has no network opcodes.  Like the file-backed editor, a program STOREs a command
# number into a control cell; the host does the work (DNS, whois, ICMP, TCP/TLS) and prints
# its output straight to stdout.  Strings are passed as the address of a NUL-terminated
# string in VM memory.  Everything has a timeout; nothing runs unless --net is given.
# ---------------------------------------------------------------------------------------
NET_MAGIC = 20037
CTL_NCMD, CTL_NA1, CTL_NA2, CTL_NA3, CTL_NA4, CTL_NSTATUS, CTL_NR1, CTL_NR2, CTL_NMAGIC = range(4112, 4121)
IANA_WHOIS = "whois.iana.org"
WHOIS_ALLOW_PRIVATE = False        # automatic whois referrals may only point at public hosts
MAX_PRINT = 256 * 1024
QTYPES = {"A": 1, "NS": 2, "CNAME": 5, "SOA": 6, "PTR": 12, "MX": 15, "TXT": 16, "AAAA": 28, "SRV": 33,
          "NAPTR": 35, "DS": 43, "DNSKEY": 48, "TLSA": 52, "CAA": 257, "ANY": 255}
QNAMES = {v: k for k, v in QTYPES.items()}
RCODES = ["NOERROR", "FORMERR", "SERVFAIL", "NXDOMAIN", "NOTIMP", "REFUSED", "YXDOMAIN", "YXRRSET",
          "NXRRSET", "NOTAUTH", "NOTZONE"]

def _esc_bytes(b):
    return "".join(chr(c) if 32 <= c < 127 and chr(c) not in '"\\' else "\\%03d" % c for c in b)

def _dns_encode_name(name):
    name = name.rstrip(".")
    if not name: return b"\0"
    out = b""
    for lab in name.split("."):
        if not lab: raise ValueError("empty label in name")
        try: b = lab.encode("ascii")
        except UnicodeEncodeError: b = lab.encode("idna")
        if len(b) > 63: raise ValueError("label longer than 63 bytes")
        out += bytes([len(b)]) + b
    if len(out) + 1 > 255: raise ValueError("name longer than 253 bytes")
    return out + b"\0"

def _dns_decode_name(msg, off):
    labels, end, hops = [], None, 0
    while True:
        if off >= len(msg): raise ValueError("truncated name")
        n = msg[off]
        if n == 0:
            off += 1; break
        if n & 0xC0 == 0xC0:
            if off + 1 >= len(msg): raise ValueError("truncated pointer")
            ptr = ((n & 0x3F) << 8) | msg[off + 1]
            if end is None: end = off + 2
            hops += 1
            if hops > 64 or ptr >= len(msg): raise ValueError("bad compression pointer")
            off = ptr; continue
        if n & 0xC0: raise ValueError("bad label type")
        off += 1
        if off + n > len(msg): raise ValueError("truncated label")
        labels.append(msg[off:off + n]); off += n
        if sum(len(l) + 1 for l in labels) > 255: raise ValueError("name too long")
    if end is None: end = off
    return ".".join(_esc_bytes(l).replace(".", "\\.") for l in labels) + ".", end

def _dns_rdata(msg, rtype, off, rdlen):
    end = off + rdlen; raw = msg[off:end]
    if rtype == 1 and rdlen == 4: return socket.inet_ntoa(raw)
    if rtype == 28 and rdlen == 16: return str(ipaddress.IPv6Address(raw))
    if rtype in (2, 5, 12): return _dns_decode_name(msg, off)[0]
    if rtype == 15 and rdlen >= 3:
        return "%d %s" % (struct.unpack(">H", raw[:2])[0], _dns_decode_name(msg, off + 2)[0])
    if rtype == 16:
        parts, i = [], 0
        while i < rdlen:
            n = raw[i]; parts.append('"' + _esc_bytes(raw[i + 1:i + 1 + n]) + '"'); i += 1 + n
        return " ".join(parts)
    if rtype == 6:
        m, o = _dns_decode_name(msg, off); r, o = _dns_decode_name(msg, o)
        return "%s %s %d %d %d %d %d" % ((m, r) + struct.unpack(">IIIII", msg[o:o + 20]))
    if rtype == 33 and rdlen >= 7:
        p, w, port = struct.unpack(">HHH", raw[:6]); return "%d %d %d %s" % (p, w, port, _dns_decode_name(msg, off + 6)[0])
    if rtype == 257 and rdlen >= 2:
        fl, tl = raw[0], raw[1]; return '%d %s "%s"' % (fl, _esc_bytes(raw[2:2 + tl]), _esc_bytes(raw[2 + tl:]))
    return "\\# %d %s" % (rdlen, raw.hex())

def _icmp_checksum(b):
    if len(b) % 2: b += b"\0"
    s = sum(struct.unpack(">%dH" % (len(b) // 2), b))
    s = (s >> 16) + (s & 0xFFFF); s += s >> 16
    return (~s) & 0xFFFF

class NetKit:
    def __init__(self, out, mem):
        self.out = out; self.mem = mem; self.sock = None

    # ------------------------------------------------------------ helpers
    def w(self, s):
        self.out.write(s); self.out.flush()
    @staticmethod
    def clean(text):
        """remote text is untrusted: no terminal escape sequences or other control characters"""
        return re.sub(r"[\x00-\x08\x0b-\x1f\x7f-\x9f]", ".", text.replace("\r", ""))
    def cs(self, addr, maxlen=255):
        if not (0 <= addr < len(self.mem)): return ""
        chars = []
        for i in range(maxlen):
            if addr + i >= len(self.mem): break
            v = self.mem[addr + i]
            if v == 0: break
            chars.append(chr(v) if 0 < v < 0x110000 else "?")
        return "".join(c for c in "".join(chars) if c.isprintable()).strip()
    @staticmethod
    def okhost(h):
        h = h.strip("[]")
        return 1 <= len(h) <= 253 and all(c.isalnum() or c in "._:%-" for c in h)
    @staticmethod
    def hostnorm(h):
        return h.strip("[]")
    def close(self):
        if self.sock is not None:
            try: self.sock.close()
            except OSError: pass
            self.sock = None

    def command(self, cmd, a1, a2, a3, a4):
        """returns (status, r1, r2)"""
        try:
            if cmd == 1: return self.c_resolve(self.cs(a1))
            if cmd == 2: return self.c_dig(self.cs(a1), self.cs(a2), self.cs(a3))
            if cmd == 3: return self.c_whois(self.cs(a1), self.cs(a2))
            if cmd == 4: return self.c_ping(self.cs(a1), a2, a3)
            if cmd == 5: return self.c_connect(self.cs(a1), a2, a3, a4)
            if cmd == 6: return self.c_send(a1, a2)
            if cmd == 7: return self.c_recv(a1, a2)
            if cmd == 8: self.close(); return (1, 0, 0)
            if cmd == 9: return self.c_listen(a1, a2, self.cs(a3))
            if cmd == 10: return self.c_http(self.cs(a1))
            if cmd == 11: return self.c_portcheck(self.cs(a1), a2, a3)
            if cmd == 12: return self.c_rdns(self.cs(a1), self.cs(a2))
            self.w("net: unknown command %s\n" % cmd); return (0, 0, 0)
        except KeyboardInterrupt:
            self.w("\n^C interrupted\n"); return (0, 0, 0)
        except Exception as e:
            self.w("error: %s: %s\n" % (type(e).__name__, e)); return (0, 0, 0)

    # ------------------------------------------------------------ resolve
    def c_resolve(self, name):
        if not self.okhost(name): self.w("error: that is not a valid host name\n"); return (0, 0, 0)
        name = self.hostnorm(name)
        try: infos = socket.getaddrinfo(name, None, 0, socket.SOCK_STREAM)
        except socket.gaierror as e:
            self.w("%s: cannot resolve (%s)\n" % (name, e.strerror or e)); return (0, 0, 0)
        seen = []
        for fam, _, _, _, sa in infos:
            ip = sa[0].split("%")[0]
            if ip not in seen: seen.append(ip)
        self.w("%s:\n" % name)
        for ip in seen: self.w("  %s  (%s)\n" % (ip, "IPv6" if ":" in ip else "IPv4"))
        return (1, len(seen), 0)

    # ------------------------------------------------------------ dig
    def _resolver(self, server):
        port = 53
        if server:
            server = server.lstrip("@")
            if "#" in server: server, p = server.rsplit("#", 1); port = int(p)
            elif server.count(":") == 1 and not server.startswith("["): server, p = server.rsplit(":", 1); port = int(p)
            server = server.strip("[]")
            if not 1 <= port <= 65535: raise ValueError("bad server port")
        else:
            server = None
            try:
                for line in open("/etc/resolv.conf"):
                    parts = line.split()
                    if len(parts) >= 2 and parts[0] == "nameserver": server = parts[1].split("%")[0]; break
            except OSError: pass
            server = server or "1.1.1.1"
        if not self.okhost(server): raise ValueError("bad server name")
        info = socket.getaddrinfo(server, port, 0, socket.SOCK_DGRAM)[0]
        return info[0], info[4]

    def _dns_udp(self, af, addr, query, qid, timeout=3.0, tries=2):
        for _ in range(tries):
            s = socket.socket(af, socket.SOCK_DGRAM)
            deadline = time.time() + timeout
            try:
                s.sendto(query, addr)
                while True:
                    left = deadline - time.time()
                    if left <= 0: break
                    s.settimeout(left)
                    try: data, src = s.recvfrom(65535)
                    except (socket.timeout, TimeoutError): break
                    if src[0].split("%")[0] != addr[0].split("%")[0]: continue
                    if len(data) >= 12 and struct.unpack(">H", data[:2])[0] == qid: return data
            finally: s.close()
        raise TimeoutError("no reply from %s (timed out)" % addr[0])

    def _dns_tcp(self, af, addr, query, qid, timeout=5.0):
        s = socket.socket(af, socket.SOCK_STREAM); s.settimeout(timeout)
        try:
            s.connect(addr); s.sendall(struct.pack(">H", len(query)) + query)
            def rd(n):
                b = b""
                while len(b) < n:
                    c = s.recv(n - len(b))
                    if not c: raise ConnectionError("connection closed by server")
                    b += c
                return b
            ln = struct.unpack(">H", rd(2))[0]; data = rd(ln)
            if len(data) < 12 or struct.unpack(">H", data[:2])[0] != qid: raise ValueError("reply id mismatch")
            return data
        finally: s.close()

    def c_dig(self, name, qtype, server):
        if qtype.startswith("@") and not server: server, qtype = qtype, ""
        qtype = (qtype or "A").upper()
        if qtype not in QTYPES: self.w("error: unknown record type %r (try %s)\n" % (qtype, " ".join(QTYPES))); return (0, 0, 0)
        if not name or " " in name: self.w("error: a name is required\n"); return (0, 0, 0)
        if qtype == "PTR":
            try: name = ipaddress.ip_address(name.strip("[]")).reverse_pointer + "."
            except ValueError: pass
        try:
            qname = _dns_encode_name(name); af, addr = self._resolver(server)
        except (ValueError, socket.gaierror) as e:
            self.w("error: %s\n" % e); return (0, 0, 0)
        qid = random.randrange(1, 65536)
        opt = b"\0" + struct.pack(">HHIH", 41, 1232, 0, 0)          # EDNS0, 1232-byte UDP payload
        query = struct.pack(">HHHHHH", qid, 0x0100, 1, 0, 0, 1) + qname + struct.pack(">HH", QTYPES[qtype], 1) + opt
        t0 = time.time(); via = "UDP"
        try:
            data = self._dns_udp(af, addr, query, qid)
            if data[2] & 0x02:
                self.w(";; Truncated, retrying in TCP mode.\n"); data = self._dns_tcp(af, addr, query, qid); via = "TCP"
        except (TimeoutError, OSError, ValueError) as e:
            self.w(";; connection error: %s\n" % e); return (0, 0, 0)
        ms = int((time.time() - t0) * 1000)
        try: self._dig_print(data, name, qtype, server, addr, ms, via)
        except (ValueError, struct.error, IndexError):
            self.w(";; error: malformed DNS reply (%d bytes)\n" % len(data)); return (0, 0, 0)
        return (1, data[3] & 0x0F, 0)

    def c_rdns(self, ip, server):
        try: ipaddress.ip_address(ip.strip("[]"))
        except ValueError: self.w("error: %r is not an IP address (use dig for names)\n" % ip); return (0, 0, 0)
        return self.c_dig(ip, "PTR", server)

    def _dig_print(self, msg, name, qtype, server, addr, ms, via):
        qid, flags, qd, an, ns, ar = struct.unpack(">HHHHHH", msg[:12])
        off, qs = 12, []
        for _ in range(qd):
            qn, off = _dns_decode_name(msg, off); qt, qc = struct.unpack(">HH", msg[off:off + 4]); off += 4; qs.append((qn, qt, qc))
        secs = []
        for count in (an, ns, ar):
            recs = []
            for _ in range(count):
                rn, off = _dns_decode_name(msg, off)
                rt, rc, ttl, rl = struct.unpack(">HHIH", msg[off:off + 10]); off += 10
                if off + rl > len(msg): raise ValueError("truncated record")
                recs.append((rn, ttl, rc, rt, _dns_rdata(msg, rt, off, rl))); off += rl
            secs.append(recs)
        rcode = flags & 0x0F
        fl = " ".join(n for bit, n in ((0x8000, "qr"), (0x0400, "aa"), (0x0200, "tc"), (0x0100, "rd"), (0x0080, "ra"), (0x0020, "ad"), (0x0010, "cd")) if flags & bit)
        self.w("; <<>> GATC dig <<>> %s %s%s\n" % (name, qtype, (" @" + server.lstrip("@")) if server else ""))
        self.w(";; ->>HEADER<<- opcode: QUERY, status: %s, id: %d\n" % (RCODES[rcode] if rcode < len(RCODES) else "RCODE%d" % rcode, qid))
        self.w(";; flags: %s; QUERY: %d, ANSWER: %d, AUTHORITY: %d, ADDITIONAL: %d\n" % (fl, qd, an, ns, ar))
        cname = lambda c: "IN" if c == 1 else "CLASS%d" % c
        tname = lambda t: QNAMES.get(t, "TYPE%d" % t)
        self.w("\n;; QUESTION SECTION:\n")
        for qn, qt, qc in qs: self.w(";%s\t\t%s\t%s\n" % (qn, cname(qc), tname(qt)))
        for title, recs in zip(("ANSWER", "AUTHORITY", "ADDITIONAL"), secs):
            recs = [r for r in recs if r[3] != 41]
            if recs:
                self.w("\n;; %s SECTION:\n" % title)
                for rn, ttl, rc, rt, rd in recs: self.w("%s\t\t%d\t%s\t%s\t%s\n" % (rn, ttl, cname(rc), tname(rt), rd))
        self.w("\n;; Query time: %d msec\n;; SERVER: %s#%d (%s)\n;; MSG SIZE  rcvd: %d\n" % (ms, addr[0].split("%")[0], addr[1], via, len(msg)))

    # ------------------------------------------------------------ whois
    @staticmethod
    def _split_server(s, default_port=43):
        if s.count(":") == 1 and not s.startswith("["):
            h, p = s.rsplit(":", 1)
            if p.isdigit(): return h, int(p)
        return s.strip("[]"), default_port
    def _public(self, host):
        if WHOIS_ALLOW_PRIVATE: return True
        try:
            for info in socket.getaddrinfo(host, 43, 0, socket.SOCK_STREAM):
                if not ipaddress.ip_address(info[4][0].split("%")[0]).is_global: return False
            return True
        except (socket.gaierror, ValueError): return False
    def _whois_fetch(self, server, query, timeout=10.0):
        host, port = self._split_server(server)
        s = socket.create_connection((host, port), timeout)
        try:
            s.sendall((query + "\r\n").encode("utf-8", "replace"))
            data, deadline = b"", time.time() + 20
            while len(data) < 512 * 1024:
                left = deadline - time.time()
                if left <= 0: break
                s.settimeout(left)
                try: chunk = s.recv(65536)
                except (socket.timeout, TimeoutError): break
                if not chunk: break
                data += chunk
        finally: s.close()
        return data.decode("utf-8", "replace")
    @staticmethod
    def _whois_fmt(server, q):
        sv = server.lower().split(":")[0]
        if sv.endswith("verisign-grs.com"): return "domain " + q
        if sv == "whois.arin.net":
            try: ipaddress.ip_address(q); return "n + " + q
            except ValueError: pass
        return q
    @staticmethod
    def _refer(text, keys):
        for key in keys:
            m = re.search(r"^\s*%s\s*:\s*(\S+)" % key, text, re.I | re.M)
            if m:
                h = m.group(1).strip().lower().rstrip(".")
                h = re.sub(r"^[a-z]+://", "", h).split("/")[0]
                if re.fullmatch(r"[a-z0-9.-]{3,253}", h): return h
        return None
    def c_whois(self, query, server):
        if not query or any(c.isspace() for c in query): self.w("error: give one domain, IP address or AS number\n"); return (0, 0, 0)
        out, path = "", []
        try:
            if server:
                path = [server]; out = self._whois_fetch(server, self._whois_fmt(server, query))
            else:
                path = [IANA_WHOIS]; out = self._whois_fetch(IANA_WHOIS, query)
                nxt = self._refer(out, ("refer", "whois"))
                for hop in range(2):
                    if not nxt or nxt in path: break
                    if not self._public(nxt):
                        self.w("# refusing to follow a referral to a non-public host (%s)\n" % nxt); break
                    got = self._whois_fetch(nxt, self._whois_fmt(nxt, query)); path.append(nxt)
                    if got.strip(): out = got
                    nxt = self._refer(got, ("Registrar WHOIS Server",)) if hop == 0 else None
        except (OSError, socket.gaierror, ValueError) as e:
            self.w("whois: %s: %s\n" % (path[-1] if path else "?", e)); return (0, 0, 0)
        lines = self.clean(out).splitlines()
        self.w("# whois: %s\n" % " -> ".join(path))
        for ln in lines[:300]: self.w(ln + "\n")
        if len(lines) > 300: self.w("# ... %d more lines not shown\n" % (len(lines) - 300))
        if not out.strip(): self.w("# (empty reply)\n"); return (0, 0, 0)
        return (1, len(lines), 0)

    # ------------------------------------------------------------ ping
    def _icmp_socket(self, fam):
        proto = socket.IPPROTO_ICMP if fam == socket.AF_INET else socket.IPPROTO_ICMPV6
        for typ, kind in ((socket.SOCK_DGRAM, "dgram"), (socket.SOCK_RAW, "raw")):
            try: return socket.socket(fam, typ, proto), kind
            except (OSError, AttributeError): continue
        return None, None
    def c_ping(self, host, count, tmo_ms):
        if not self.okhost(host): self.w("error: that is not a valid host name\n"); return (0, 0, 0)
        host = self.hostnorm(host)
        count = max(1, min(count if count > 0 else 4, 100)); tmo = max(200, min(tmo_ms if tmo_ms > 0 else 2000, 10000)) / 1000.0
        try: infos = socket.getaddrinfo(host, None, 0, socket.SOCK_STREAM)
        except socket.gaierror as e: self.w("ping: %s: cannot resolve (%s)\n" % (host, e.strerror or e)); return (0, 0, 0)
        infos.sort(key=lambda i: i[0] != socket.AF_INET)
        fam, _, _, _, sa = infos[0]; ip = sa[0].split("%")[0]; v4 = fam == socket.AF_INET
        sock, kind = self._icmp_socket(fam)
        if sock is None: return self._tcp_ping(host, ip, count, tmo)
        self.w("PING %s (%s): 56 data bytes\n" % (host, ip))
        ident, rtts, sent = os.getpid() & 0xFFFF, [], 0
        want = 0 if v4 else 129
        try:
            for seq in range(1, count + 1):
                payload = struct.pack(">d", time.time()) + bytes(range(48))
                hdr = struct.pack(">BBHHH", 8 if v4 else 128, 0, 0, ident, seq)
                pkt = hdr + payload
                if v4: pkt = struct.pack(">BBHHH", 8, 0, _icmp_checksum(hdr + payload), ident, seq) + payload
                t0 = time.time(); sent += 1
                try: sock.sendto(pkt, sa)
                except OSError as e: self.w("ping: send failed: %s\n" % e); break
                deadline = t0 + tmo; got = False
                while True:
                    left = deadline - time.time()
                    if left <= 0: break
                    sock.settimeout(left)
                    try: data, src = sock.recvfrom(2048)
                    except (socket.timeout, TimeoutError): break
                    ttl = None
                    if v4 and data and data[0] >> 4 == 4:
                        ihl = (data[0] & 0x0F) * 4; ttl = data[8]; data = data[ihl:]
                    if len(data) < 8: continue
                    typ, _, _, rid, rseq = struct.unpack(">BBHHH", data[:8])
                    if typ != want or rseq != seq or (kind == "raw" and rid != ident): continue
                    rtt = (time.time() - t0) * 1000; rtts.append(rtt); got = True
                    self.w("%d bytes from %s: icmp_seq=%d%s time=%.2f ms\n" % (len(data), ip, seq, (" ttl=%d" % ttl) if ttl else "", rtt))
                    break
                if not got: self.w("Request timeout for icmp_seq %d\n" % seq)
                if seq < count:
                    rest = 1.0 - (time.time() - t0)
                    if rest > 0: time.sleep(rest)
        finally: sock.close()
        return self._ping_summary(host, sent, rtts)
    def _ping_summary(self, host, sent, rtts):
        loss = 100.0 * (sent - len(rtts)) / sent if sent else 100.0
        self.w("\n--- %s ping statistics ---\n%d packets transmitted, %d received, %.0f%% packet loss\n" % (host, sent, len(rtts), loss))
        if rtts:
            avg = sum(rtts) / len(rtts); md = (sum((r - avg) ** 2 for r in rtts) / len(rtts)) ** 0.5
            self.w("rtt min/avg/max/mdev = %.3f/%.3f/%.3f/%.3f ms\n" % (min(rtts), avg, max(rtts), md))
        return (1 if rtts else 0, len(rtts), sent)
    def _tcp_ping(self, host, ip, count, tmo):
        self.w("PING %s (%s): ICMP is not permitted here, using TCP connects instead\n" % (host, ip))
        port, rtts, sent = 443, [], 0
        for seq in range(1, count + 1):
            sent += 1; t0 = time.time(); res = None
            for attempt in range(2):
                try:
                    s = socket.create_connection((ip, port), tmo); s.close(); res = "open"
                except ConnectionRefusedError: res = "refused"
                except OSError: res = None
                if res or port == 80: break
                port = 80
            if res:
                rtt = (time.time() - t0) * 1000; rtts.append(rtt)
                self.w("tcp %s:%d seq=%d %s time=%.2f ms\n" % (ip, port, seq, res, rtt))
            else: self.w("tcp %s:%d seq=%d no response\n" % (ip, port, seq))
            if seq < count: time.sleep(max(0, 1.0 - (time.time() - t0)))
        return self._ping_summary(host, sent, rtts)

    # ------------------------------------------------------------ tcp client / server (nc)
    def _tls_info(self, ss):
        info = "%s %s" % (ss.version(), (ss.cipher() or ("?",))[0])
        try:
            cert = ss.getpeercert(); subj = dict(x[0] for x in cert.get("subject", ())); iss = dict(x[0] for x in cert.get("issuer", ()))
            left = (ssl.cert_time_to_seconds(cert["notAfter"]) - time.time()) / 86400.0
            return info, "cert: %s  issuer: %s  expires in %d days" % (subj.get("commonName", "?"), iss.get("organizationName") or iss.get("commonName", "?"), left)
        except Exception: return info, None
    def _open(self, host, port, timeout, tls):
        """returns (socket, connect_ms, ip, tls_info) or raises"""
        last = None
        for fam, typ, proto, _, sa in socket.getaddrinfo(host, port, 0, socket.SOCK_STREAM):
            s = socket.socket(fam, typ, proto); s.settimeout(timeout); t0 = time.time()
            try:
                s.connect(sa); ms = (time.time() - t0) * 1000
                if tls:
                    s = ssl.create_default_context().wrap_socket(s, server_hostname=host)
                return s, ms, sa[0].split("%")[0]
            except (OSError, ssl.SSLError) as e:
                last = e
                try: s.close()
                except OSError: pass
        raise last
    def _why(self, e, host):
        if isinstance(e, socket.gaierror): return "cannot resolve %s (%s)" % (host, e.strerror or e)
        if isinstance(e, ssl.SSLCertVerificationError): return "TLS certificate check failed: %s" % (e.verify_message or e)
        if isinstance(e, ssl.SSLError): return "TLS error: %s" % e
        if isinstance(e, (socket.timeout, TimeoutError)): return "timed out"
        if isinstance(e, ConnectionRefusedError): return "connection refused"
        return e.strerror or str(e)
    def c_connect(self, host, port, tmo, tls):
        self.close()
        if not self.okhost(host): self.w("error: that is not a valid host name\n"); return (0, 0, 0)
        if not 1 <= port <= 65535: self.w("error: port must be 1-65535\n"); return (0, 0, 0)
        host = self.hostnorm(host); tmo = max(1, min(tmo or 5, 60))
        try: s, ms, ip = self._open(host, port, tmo, bool(tls))
        except (OSError, ssl.SSLError) as e: self.w("connect to %s port %d failed: %s\n" % (host, port, self._why(e, host))); return (0, 0, 0)
        self.sock = s
        self.w("connected to %s (%s) port %d in %.0f ms\n" % (host, ip, port, ms))
        if tls:
            info, cert = self._tls_info(s); self.w("tls: %s\n" % info)
            if cert: self.w("%s\n" % cert)
        return (1, int(ms), 0)
    def c_send(self, addr, n):
        if self.sock is None: self.w("error: not connected\n"); return (0, 0, 0)
        if not (0 <= n <= 4096 and 0 <= addr and addr + n <= len(self.mem)): self.w("error: bad send buffer\n"); return (0, 0, 0)
        try: self.sock.settimeout(10); self.sock.sendall(bytes(self.mem[addr + i] & 0xFF for i in range(n)))
        except (OSError, ssl.SSLError) as e:
            self.w("send failed: %s\n" % self._why(e, "")); self.close(); return (0, 0, 0)
        return (1, n, 0)
    def c_recv(self, idle_ms, total_ms):
        if self.sock is None: self.w("error: not connected\n"); return (0, 0, 1)
        idle = max(0.05, min(idle_ms or 1500, 60000) / 1000.0); total = max(idle, min(total_ms or 5000, 300000) / 1000.0)
        deadline, data, closed = time.time() + total, b"", 0
        while len(data) < MAX_PRINT:
            left = min(idle, deadline - time.time())
            if left <= 0: break
            self.sock.settimeout(left)
            try: chunk = self.sock.recv(65536)
            except (socket.timeout, TimeoutError): break
            except ssl.SSLWantReadError: continue
            except (OSError, ssl.SSLError): closed = 1; break
            if not chunk: closed = 1; break
            data += chunk
        if data:
            text = self.clean(data.decode("utf-8", "replace"))
            self.w(text if text.endswith("\n") else text + "\n")
        if closed: self.close()
        return (1 if data else 0, len(data), closed)
    def c_listen(self, port, tmo, bind):
        self.close()
        bind = bind or "127.0.0.1"
        if not self.okhost(bind): self.w("error: bad bind address\n"); return (0, 0, 0)
        if not 1 <= port <= 65535: self.w("error: port must be 1-65535\n"); return (0, 0, 0)
        tmo = max(1, min(tmo or 60, 3600)); bind = self.hostnorm(bind)
        try:
            fam = socket.getaddrinfo(bind, port, 0, socket.SOCK_STREAM, 0, socket.AI_PASSIVE)[0]
            srv = socket.socket(fam[0], socket.SOCK_STREAM); srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            srv.bind(fam[4]); srv.listen(1); srv.settimeout(tmo)
        except (OSError, socket.gaierror) as e:
            self.w("listen on %s port %d failed: %s\n" % (bind, port, self._why(e, bind))); return (0, 0, 0)
        self.w("listening on %s port %d (waits %d s, Ctrl-C cancels)\n" % (bind, port, tmo))
        try: conn, peer = srv.accept()
        except (socket.timeout, TimeoutError): self.w("no connection received (timed out)\n"); return (0, 0, 0)
        finally: srv.close()
        self.sock = conn; self.w("connection from %s port %d\n" % (peer[0].split("%")[0], peer[1]))
        return (1, 0, 0)
    def c_portcheck(self, host, port, tmo):
        if not self.okhost(host): self.w("error: that is not a valid host name\n"); return (0, 0, 0)
        if not 1 <= port <= 65535: self.w("port %d: bad port number\n" % port); return (0, 0, 0)
        host = self.hostnorm(host); tmo = max(1, min(tmo or 3, 30))
        try: s, ms, ip = self._open(host, port, tmo, False); s.close()
        except socket.gaierror as e: self.w("%s: cannot resolve (%s)\n" % (host, e.strerror or e)); return (0, 0, 0)
        except (socket.timeout, TimeoutError): self.w("%s %d/tcp: no response (timed out - filtered or host down)\n" % (host, port)); return (0, 0, 0)
        except ConnectionRefusedError: self.w("%s %d/tcp: closed (connection refused)\n" % (host, port)); return (0, 0, 0)
        except OSError as e: self.w("%s %d/tcp: %s\n" % (host, port, e.strerror or e)); return (0, 0, 0)
        self.w("%s %d/tcp: open (%s, %.0f ms)\n" % (host, port, ip, ms)); return (1, int(ms), 0)

    # ------------------------------------------------------------ http HEAD
    def c_http(self, url):
        from urllib.parse import urlsplit, quote
        u = url if "://" in url else "http://" + url
        try: sp = urlsplit(u); port = sp.port
        except ValueError: self.w("error: bad URL\n"); return (0, 0, 0)
        scheme = sp.scheme.lower()
        if scheme not in ("http", "https") or not sp.hostname or not self.okhost(sp.hostname):
            self.w("error: use http://host/path or https://host/path\n"); return (0, 0, 0)
        tls = scheme == "https"; port = port or (443 if tls else 80)
        path = quote(sp.path or "/", safe="/%:@!$&'()*+,;=~-._") + (("?" + quote(sp.query, safe="=&%:@!$'()*+,;/?~-._")) if sp.query else "")
        hosthdr = sp.hostname + ((":%d" % port) if port != (443 if tls else 80) else "")
        self.w("HEAD %s://%s%s\n" % (scheme, hosthdr, path))
        t0 = time.time()
        try: s, ms, ip = self._open(sp.hostname, port, 8, tls)
        except (OSError, ssl.SSLError) as e: self.w("failed: %s\n" % self._why(e, sp.hostname)); return (0, 0, 0)
        try:
            self.w("connected to %s port %d in %.0f ms\n" % (ip, port, ms))
            if tls:
                info, cert = self._tls_info(s); self.w("tls: %s\n" % info)
                if cert: self.w("%s\n" % cert)
            s.settimeout(8)
            s.sendall(("HEAD %s HTTP/1.1\r\nHost: %s\r\nUser-Agent: GATC-netkit/1.0\r\nAccept: */*\r\nConnection: close\r\n\r\n" % (path, hosthdr)).encode("ascii", "replace"))
            buf, first = b"", None
            while b"\r\n\r\n" not in buf and len(buf) < 65536:
                try: chunk = s.recv(8192)
                except (socket.timeout, TimeoutError): break
                if not chunk: break
                if first is None: first = (time.time() - t0) * 1000
                buf += chunk
        except (OSError, ssl.SSLError) as e: self.w("failed: %s\n" % self._why(e, sp.hostname)); return (0, 0, 0)
        finally: s.close()
        head = buf.split(b"\r\n\r\n")[0].decode("latin-1")
        if not head.strip(): self.w("no response\n"); return (0, 0, 0)
        self.w("\n%s\n" % self.clean(head))
        if first is not None: self.w("\ntime to first byte: %.0f ms\n" % first)
        m = re.match(r"HTTP/\S+\s+(\d+)", head)
        return (1, int(m.group(1)) if m else 0, 0)

# ---------------------------------------------------------------------------------------
# filesystem bridge for sak.aa  (opt-in: run with --fs)
#
# The VM has no file-system opcodes.  Exactly like the --net bridge, a program STOREs a
# command number into a control cell and the host performs ONE small syscall-like step
# (stat, open, read, write, unlink, rename, ...).  Everything else - parsing, globbing,
# recursion, formatting, grep/wc/tail/xxd/crc32 - is done by the .aa program itself.
# Paths are the address of a NUL-terminated byte string in VM memory (one byte per cell).
# Nothing here is reachable unless --fs is given.  --fs-root DIR confines it to DIR.
# Control cells 4128..4138:
#   4128 command  4129-4132 arguments  4133 status (1 ok, 0 failed)
#   4134-4137 results (on failure 4134 = errno)  4138 magic (FS_MAGIC on, 1 = supported but off)
# Commands:
#   1 stat(path,follow)->size,type,mode,mtime   type: 1 file 2 dir 3 symlink 4 other
#   2 opendir(path)->handle,count   3 dirent(h,i,dest,max)->len,type,mode,size   4 closedir(h)
#   5 mkdir(path,parents)  6 rmdir(path)  7 unlink(path)  8 rename(src,dst)
#   9 open(path,mode)->fd   mode: 0 read 1 write+truncate 2 append 3 create (no truncate)
#   10 read(fd,addr,n)->count  11 write(fd,addr,n)->count  (fd 1 = stdout, 2 = stderr)
#   12 close(fd)  13 seek(fd,off,whence,negative)->pos  14 chdir(path)  15 getcwd(dest,max)->len
#   16 utime(path)  17 chmod(path,mode)  18 time()->epoch,utc offset  19 strerror(errno,dest,max)
#   20 redirect stdout(path or 0 to restore, append)  21 isatty()->stdout tty,stdin tty,columns
#   22 realpath(path,dest,max)  23 getenv(name,dest,max)  24 samefile(a,b)->same,b inside a
#   25 errmode(1 = PUTC/PUTN/raw writes go to stderr, 0 = back)  26 statvfs(path)->total,free,avail
#   27 link(target,linkpath,hard)  28 readlink(path,dest,max)->len
# ---------------------------------------------------------------------------------------
FS_MAGIC = 18003
CTL_FCMD, CTL_FA1, CTL_FA2, CTL_FA3, CTL_FA4, CTL_FSTATUS, CTL_FR1, CTL_FR2, CTL_FR3, CTL_FR4, CTL_FMAGIC = range(4128, 4139)
FS_IO_MAX = 4096
FS_MAX_HANDLES = 64

class FSKit:
    def __init__(self, vm, root=None):
        import io, stat as _stat
        self.io, self.st = io, _stat
        self.vm = vm; self.mem = vm.mem
        self.fds = {}; self.dirs = {}; self.next_id = 3
        self.term_out = vm.stdout           # the real stdout text stream
        self.redir_bin = None; self.redir_out = None
        self.err_mode = False
        self.root = os.path.realpath(root) if root else None
        if self.root:
            if not os.path.isdir(self.root): raise GATCError(f"--fs-root {root}: not a directory")
            os.chdir(self.root)

    # ------------------------------------------------------------ helpers
    def cs(self, addr, maxlen=FS_IO_MAX):
        """NUL-terminated byte string at addr -> str path (round-trips any byte via fsdecode)"""
        if not (0 <= addr < len(self.mem)): raise OSError(14, "bad address")
        b = bytearray()
        for i in range(maxlen):
            if addr + i >= len(self.mem): break
            v = self.mem[addr + i]
            if v == 0: break
            b.append(v & 0xFF)
        s = os.fsdecode(bytes(b))
        if s == "~" or s.startswith("~/"):
            home = self.root or os.path.expanduser("~")
            s = home + s[1:]
        return s

    def put(self, dest, maxlen, data):
        """write bytes + NUL into VM memory, truncated to fit; returns the length written"""
        if maxlen <= 0 or not (0 <= dest < len(self.mem)): return 0
        maxlen = min(maxlen, len(self.mem) - dest)
        data = data[:maxlen - 1]
        for i, c in enumerate(data): self.mem[dest + i] = c
        self.mem[dest + len(data)] = 0
        return len(data)

    def jail(self, path, follow=True):
        """with --fs-root, refuse paths that resolve outside the root (best effort)"""
        if not self.root: return path
        full = os.path.abspath(path)
        if follow: real = os.path.realpath(full)
        else: real = os.path.join(os.path.realpath(os.path.dirname(full)), os.path.basename(full))
        if real != self.root and not real.startswith(self.root.rstrip(os.sep) + os.sep):
            raise OSError(13, "outside --fs-root")
        return path

    def ftype(self, mode):
        st = self.st
        if st.S_ISLNK(mode): return 3
        if st.S_ISDIR(mode): return 2
        if st.S_ISREG(mode): return 1
        return 4

    def new_id(self, table):
        if len(self.fds) + len(self.dirs) >= FS_MAX_HANDLES: raise OSError(24, "too many open handles")
        h = self.next_id; self.next_id += 1
        return h

    def cur_out(self):
        return self.redir_out if self.redir_out is not None else self.term_out

    def raw_out(self, data, fd):
        """byte-exact write to stdout (fd 1) or stderr (fd 2), in order with PUTC/PUTN output"""
        if fd == 2 or self.err_mode:
            sys.stderr.flush(); target = getattr(sys.stderr, "buffer", None)
            if target is None: sys.stderr.write(data.decode("latin-1"))
            else: target.write(data); target.flush()
            return
        if self.redir_out is not None:
            self.redir_out.flush(); self.redir_bin.write(data); self.redir_bin.flush(); return
        out = self.term_out; out.flush()
        buf = getattr(out, "buffer", None)
        if buf is None: out.write(data.decode("utf-8", "replace"))
        else: buf.write(data); buf.flush()

    def restore_stdout(self):
        if self.redir_out is not None:
            try: self.redir_out.flush(); self.redir_out.close()
            except (OSError, ValueError): pass
            self.redir_out = self.redir_bin = None
        self.vm.stdout = sys.stderr if self.err_mode else self.term_out

    def close_all(self):
        for fd in list(self.fds.values()):
            try: os.close(fd)
            except OSError: pass
        self.fds.clear(); self.dirs.clear()
        self.err_mode = False
        self.restore_stdout()

    # ------------------------------------------------------------ dispatcher
    def command(self, cmd, a1, a2, a3, a4):
        """returns (status, r1, r2, r3, r4); on failure r1 = errno"""
        try:
            f = getattr(self, "c%d" % cmd, None)
            if f is None: return (0, 22, 0, 0, 0)
            r = f(a1, a2, a3, a4)
            return (1,) + tuple(r) + (0,) * (4 - len(r))
        except OSError as e:
            return (0, e.errno or 5, 0, 0, 0)
        except (ValueError, UnicodeError, OverflowError) as e:
            return (0, 22, 0, 0, 0)

    def c1(self, a1, a2, a3, a4):     # stat
        p = self.jail(self.cs(a1), follow=bool(a2))
        s = os.stat(p) if a2 else os.lstat(p)
        return (s.st_size, self.ftype(s.st_mode), self.st.S_IMODE(s.st_mode), int(s.st_mtime))
    def c2(self, a1, a2, a3, a4):     # opendir
        p = self.jail(self.cs(a1))
        names = os.listdir(p)
        names.sort(key=lambda n: (n.lstrip(".").lower(), n))
        h = self.new_id(self.dirs); self.dirs[h] = (p, names)
        return (h, len(names))
    def c3(self, a1, a2, a3, a4):     # dirent
        p, names = self.dirs[a1] if a1 in self.dirs else (None, None)
        if names is None: raise OSError(9, "bad directory handle")
        if not 0 <= a2 < len(names): raise OSError(22, "index out of range")
        name = names[a2]
        try:
            s = os.lstat(os.path.join(p, name)); t, m, z = self.ftype(s.st_mode), self.st.S_IMODE(s.st_mode), s.st_size
        except OSError: t, m, z = 4, 0, 0
        n = self.put(a3, a4, os.fsencode(name))
        return (n, t, m, z)
    def c4(self, a1, a2, a3, a4):     # closedir
        self.dirs.pop(a1, None); return ()
    def c5(self, a1, a2, a3, a4):     # mkdir
        p = self.jail(self.cs(a1), follow=False)
        if a2: os.makedirs(p, exist_ok=True)
        else: os.mkdir(p)
        return ()
    def c6(self, a1, a2, a3, a4):     # rmdir
        p = self.jail(self.cs(a1), follow=False)
        if os.path.abspath(p) == os.path.abspath(os.sep) or (self.root and os.path.realpath(p) == self.root):
            raise OSError(16, "refusing to remove the root")
        os.rmdir(p); return ()
    def c7(self, a1, a2, a3, a4):     # unlink
        os.unlink(self.jail(self.cs(a1), follow=False)); return ()
    def c8(self, a1, a2, a3, a4):     # rename
        os.rename(self.jail(self.cs(a1), follow=False), self.jail(self.cs(a2), follow=False)); return ()
    def c9(self, a1, a2, a3, a4):     # open
        flags = {0: os.O_RDONLY, 1: os.O_WRONLY | os.O_CREAT | os.O_TRUNC,
                 2: os.O_WRONLY | os.O_CREAT | os.O_APPEND, 3: os.O_WRONLY | os.O_CREAT}.get(a2)
        if flags is None: raise OSError(22, "bad open mode")
        p = self.jail(self.cs(a1), follow=True)
        h = self.new_id(self.fds)
        self.fds[h] = os.open(p, flags | getattr(os, "O_BINARY", 0), 0o666)
        return (h,)
    def _buf(self, addr, n):
        if not (0 <= n <= FS_IO_MAX and 0 <= addr and addr + n <= len(self.mem)): raise OSError(14, "bad buffer")
    def c10(self, a1, a2, a3, a4):    # read
        self._buf(a2, a3)
        if a1 not in self.fds: raise OSError(9, "bad file handle")
        data = os.read(self.fds[a1], a3)
        for i, c in enumerate(data): self.mem[a2 + i] = c
        return (len(data),)
    def c11(self, a1, a2, a3, a4):    # write
        self._buf(a2, a3)
        data = bytes(self.mem[a2 + i] & 0xFF for i in range(a3))
        if a1 in (1, 2): self.raw_out(data, a1); return (len(data),)
        if a1 not in self.fds: raise OSError(9, "bad file handle")
        mv, done = memoryview(data), 0
        while done < len(data): done += os.write(self.fds[a1], mv[done:])
        return (done,)
    def c12(self, a1, a2, a3, a4):    # close
        fd = self.fds.pop(a1, None)
        if fd is None: raise OSError(9, "bad file handle")
        os.close(fd); return ()
    def c13(self, a1, a2, a3, a4):    # seek
        if a1 not in self.fds: raise OSError(9, "bad file handle")
        return (os.lseek(self.fds[a1], -a2 if a4 else a2, a3),)
    def c14(self, a1, a2, a3, a4):    # chdir
        os.chdir(self.jail(self.cs(a1))); return ()
    def c15(self, a1, a2, a3, a4):    # getcwd
        return (self.put(a1, a2, os.fsencode(os.getcwd())),)
    def c16(self, a1, a2, a3, a4):    # utime
        os.utime(self.jail(self.cs(a1)), None); return ()
    def c17(self, a1, a2, a3, a4):    # chmod
        os.chmod(self.jail(self.cs(a1)), a2 & 0o7777); return ()
    def c18(self, a1, a2, a3, a4):    # time
        now = int(time.time())
        return (now, time.localtime(now).tm_gmtoff or 0)
    def c19(self, a1, a2, a3, a4):    # strerror
        return (self.put(a2, a3, os.strerror(a1).encode("utf-8", "replace")),)
    def c20(self, a1, a2, a3, a4):    # redirect stdout
        self.vm.stdout.flush()
        if a1 == 0: self.restore_stdout(); return ()
        p = self.jail(self.cs(a1), follow=True)
        f = open(p, "ab" if a2 else "wb")
        if self.redir_out is not None: self.restore_stdout()
        self.redir_bin = f
        self.redir_out = self.io.TextIOWrapper(f, encoding="latin-1", newline="\n", write_through=True)
        if not self.err_mode: self.vm.stdout = self.redir_out
        return ()
    def c21(self, a1, a2, a3, a4):    # isatty
        def tty(s):
            try: return 1 if s.isatty() else 0
            except (AttributeError, ValueError): return 0
        out = 0 if self.redir_out is not None else tty(self.term_out)
        try:
            import shutil; cols = shutil.get_terminal_size((80, 24)).columns
        except (OSError, ValueError): cols = 80
        return (out, tty(sys.stdin), max(20, min(cols, 400)))
    def c22(self, a1, a2, a3, a4):    # realpath
        return (self.put(a2, a3, os.fsencode(os.path.realpath(self.jail(self.cs(a1))))),)
    def c23(self, a1, a2, a3, a4):    # getenv
        name = self.cs(a1)
        if name == "HOME" and self.root: v = self.root
        else: v = os.environ.get(name)
        if v is None: raise OSError(2, "not set")
        return (self.put(a2, a3, os.fsencode(v)),)
    def c24(self, a1, a2, a3, a4):    # samefile + inside
        a, b = self.cs(a1), self.cs(a2)
        try: same = 1 if os.path.samefile(a, b) else 0
        except OSError: same = 0
        ra, rb = os.path.realpath(a), os.path.realpath(b)
        inside = 1 if rb.startswith(ra.rstrip(os.sep) + os.sep) else 0
        return (same, inside)
    def c25(self, a1, a2, a3, a4):    # errmode
        self.vm.stdout.flush()
        self.err_mode = bool(a1)
        self.vm.stdout = sys.stderr if self.err_mode else self.cur_out()
        return ()
    def c26(self, a1, a2, a3, a4):    # statvfs
        s = os.statvfs(self.jail(self.cs(a1)))
        return (s.f_blocks * s.f_frsize, s.f_bfree * s.f_frsize, s.f_bavail * s.f_frsize)
    def c27(self, a1, a2, a3, a4):    # link / symlink
        target, link = self.cs(a1), self.jail(self.cs(a2), follow=False)
        if a3: os.link(self.jail(target), link)
        else: os.symlink(target, link)
        return ()
    def c28(self, a1, a2, a3, a4):    # readlink
        return (self.put(a2, a3, os.fsencode(os.readlink(self.jail(self.cs(a1), follow=False)))),)

class VM:
    def __init__(self, program, genome=False, stdin=None, stdout=None):
        self.program=program; self.pc=0; self.stack=[]; self.calls=[]; self.mem=[0]*MEM_SIZE
        self.stdin=stdin or sys.stdin; self.stdout=stdout or sys.stdout
        self.halted=False; self.exit_code=0; self.steps=0; self.genome=genome; self._halt_note=False
        self.file=bytearray(); self.filepath=None; self.savepath=None; self.dirty=0
        self.ctl={k:0 for k in range(CTL_A1,CTL_LOST+1)}; self.ctl.update({k:0 for k in range(CTL_NA1,CTL_NR2+1)}); self.net=None
        self.ctl.update({k:0 for k in range(CTL_FA1,CTL_FR4+1)}); self.fs=None
        self.groups=collections.deque(); self.u_total=0; self.u_suspend=0
    def load_file(self, path, addr=0):
        try:
            data = open(path, "rb").read()
        except OSError as e:
            print(f"error: cannot open {path}: {e}", file=sys.stderr)
            return
        n = min(len(data), len(self.mem) - addr)
        for i in range(n):
            self.mem[addr + i] = data[i]
        print(f"loaded {n} bytes from {path} at address {addr}", file=sys.stderr)
    # ---------------- file-backed address space ----------------
    def load_edit_file(self, path):
        if not os.path.exists(path):
            self.filepath = path; print(f"new file {path} (empty)", file=sys.stderr); return
        try:
            size = os.path.getsize(path)
            if size > MAX_FILE:
                raise GATCError(f"{path} is {size} bytes, over the 1 GiB limit ({MAX_FILE}) - not opened")
            buf = bytearray(size); pos = 0
            with open(path, "rb") as f:
                mv = memoryview(buf)
                while pos < size:
                    got = f.readinto(mv[pos:pos + CHUNK])
                    if not got: break
                    pos += got
                del mv
            if pos != size: raise GATCError(f"short read on {path}: {pos} of {size} bytes")
        except OSError as e:
            raise GATCError(f"cannot open {path}: {e}")
        self.file = buf; self.filepath = path
        print(f"opened {path}: {len(self.file)} bytes", file=sys.stderr)
    def fetch(self, addr):
        if 0 <= addr < len(self.mem): return self.mem[addr]
        if addr >= FILE_BASE:
            i = addr - FILE_BASE
            return self.file[i] if i < len(self.file) else 0
        if addr == CTL_SIZE: return len(self.file)
        if addr == CTL_DIRTY: return self.dirty
        if addr == CTL_MAGIC: return 51966
        if addr == CTL_NMAGIC: return NET_MAGIC if self.net is not None else 1   # 1 = supported but off (no --net)
        if addr == CTL_FMAGIC: return FS_MAGIC if self.fs is not None else 1     # 1 = supported but off (no --fs)
        return self.ctl.get(addr, 0)
    def store(self, addr, val):
        if 0 <= addr < len(self.mem): self.mem[addr] = val; return
        if addr >= FILE_BASE:
            i, v = addr - FILE_BASE, val & 0xFF
            n = len(self.file)
            if i < n:
                if self.file[i] != v:
                    self.u_rec(('r', i, bytes((self.file[i],)), n), 120)
                    self.file[i] = v; self.dirty = 1
            elif i < MAX_FILE and i - n <= GROW_CAP:
                self._prep(n, i - n + 1); self.file[i] = v; self.dirty = 1
            else: warn(f"file store too far past the end: {i}")
            return
        if addr == CTL_CMD: self.ctl_command(val); return
        if addr in (CTL_A1, CTL_A2, CTL_A3): self.ctl[addr] = val; return
        if addr == CTL_NCMD: self.net_command(val); return
        if addr in (CTL_NA1, CTL_NA2, CTL_NA3, CTL_NA4): self.ctl[addr] = val; return
        if addr == CTL_FCMD: self.fs_command(val); return
        if addr in (CTL_FA1, CTL_FA2, CTL_FA3, CTL_FA4): self.ctl[addr] = val; return
        warn(f"store out of bounds {addr}")
    def enable_net(self):
        self.net = NetKit(self.stdout, self.mem)
    def enable_fs(self, root=None):
        self.fs = FSKit(self, root)
    def fs_command(self, cmd):
        c = self.ctl
        if self.fs is None: c[CTL_FSTATUS] = 0; c[CTL_FR1] = 1; return     # EPERM: --fs not given
        r = self.fs.command(cmd, c[CTL_FA1], c[CTL_FA2], c[CTL_FA3], c[CTL_FA4])
        c[CTL_FSTATUS], c[CTL_FR1], c[CTL_FR2], c[CTL_FR3], c[CTL_FR4] = r
        if self.net is not None: self.net.out = self.stdout
    def net_command(self, cmd):
        c = self.ctl
        if self.net is None: c[CTL_NSTATUS] = 0; return
        st, r1, r2 = self.net.command(cmd, c[CTL_NA1], c[CTL_NA2], c[CTL_NA3], c[CTL_NA4])
        c[CTL_NSTATUS], c[CTL_NR1], c[CTL_NR2] = st, r1, r2
    # ---------------- undo history (kept in the host) ----------------
    def u_begin(self):
        self.u_suspend = 0
        if not self.groups or self.groups[-1][0]: self.groups.append([[], 0])
    def u_drop(self):
        self.groups.clear(); self.u_total = 0; self.u_suspend = 1; self.ctl[CTL_LOST] = 1
    def u_rec(self, op, cost):
        if self.u_suspend: return
        if not self.groups: self.groups.append([[], 0])
        g = self.groups[-1]; g[0].append(op); g[1] += cost; self.u_total += cost
        while self.u_total > UNDO_CAP and len(self.groups) > 1:
            self.u_total -= self.groups.popleft()[1]
        if self.u_total > UNDO_CAP: self.u_drop()
    def u_undo(self):
        while self.groups and not self.groups[-1][0]: self.groups.pop()
        if not self.groups: return 0
        ops, cost = self.groups.pop(); self.u_total -= cost
        for op in reversed(ops):
            k = op[0]
            if k == 'r':
                _, off, old, oldlen = op
                self.file[off:off + len(old)] = old
                if len(self.file) > oldlen: del self.file[oldlen:]
            elif k == 'i': del self.file[op[1]:op[1] + op[2]]
            elif k == 'd': self.file[op[1]:op[1]] = op[2]
        self.dirty = 1; self.ctl[CTL_R1] = len(ops); self.ctl[CTL_LOST] = 0
        return 1
    def _prep(self, off, ln):
        """record the old bytes of file[off:off+ln], then zero-extend the file to cover the range"""
        n = len(self.file); op = None
        if ln <= 0: return None
        if ln > UNDO_CAP // 2: self.u_drop()
        else:
            old = bytes(self.file[off:min(off + ln, n)]) if off < n else b""
            op = ('r', off, old, n); self.u_rec(op, 120 + len(old))
        while n < off + ln:
            step = min(CHUNK, off + ln - n); self.file.extend(bytes(step)); n += step
        self.dirty = 1
        return op
    def _commit(self, op, off, ln):
        """if the operation changed nothing, take its undo record back out again"""
        if op is None or not self.groups: return
        ops, cost = self.groups[-1][0], self.groups[-1][1]
        if ops and ops[-1] is op and op[3] == len(self.file) and op[2] == bytes(self.file[off:off + ln]):
            ops.pop(); self.groups[-1][1] -= 120 + len(op[2]); self.u_total -= 120 + len(op[2])
    # ---------------- bulk commands ----------------
    def ctl_command(self, cmd):
        c = self.ctl; off, cnt, a3 = c[CTL_A1], c[CTL_A2], c[CTL_A3]
        n = len(self.file); c[CTL_STATUS] = 0
        if cmd == 7: self.u_begin(); c[CTL_STATUS] = 1
        elif cmd == 6: c[CTL_STATUS] = self.u_undo()
        elif cmd == 1:
            path = self.savepath or self.filepath
            if not path: print("error: nothing to save to - use --file FILE or --save FILE", file=sys.stderr); return
            tmp = path + ".gatc-tmp"
            try:
                with open(tmp, "wb") as f: f.write(self.file)
                if os.path.exists(path):
                    import shutil; shutil.copymode(path, tmp)
                os.replace(tmp, path)
            except OSError as e:
                print(f"error: cannot write {path}: {e}", file=sys.stderr); return
            self.dirty = 0; c[CTL_STATUS] = 1
        elif cmd == 2:
            if 0 <= off <= n and 0 <= cnt <= GROW_CAP and n + cnt <= MAX_FILE:
                if cnt: self.u_rec(('i', off, cnt), 120); self.file[off:off] = bytes(cnt); self.dirty = 1
                c[CTL_R1] = cnt; c[CTL_STATUS] = 1
        elif cmd == 3:
            if 0 <= off < n and cnt >= 0:
                cnt = min(cnt, n - off); removed = bytes(self.file[off:off + cnt])
                if cnt:
                    if len(removed) > UNDO_CAP // 2: self.u_drop()
                    else: self.u_rec(('d', off, removed), 120 + len(removed))
                    del self.file[off:off + cnt]; self.dirty = 1
                c[CTL_R1] = cnt; c[CTL_STATUS] = 1
        elif cmd == 4:
            if 0 < cnt <= 64 and 0 <= a3 and a3 + cnt <= len(self.mem) and off >= 0:
                pat = bytes(v & 0xFF for v in self.mem[a3:a3 + cnt]); i = self.file.find(pat, off)
                if i >= 0: c[CTL_R1] = i; c[CTL_STATUS] = 1
        elif cmd == 5:
            if 0 <= off <= n and 0 <= cnt and off + cnt <= n:
                s = x = 0; crc = 0
                for p in range(off, off + cnt, CHUNK):
                    part = self.file[p:min(p + CHUNK, off + cnt)]
                    s += sum(part); crc = zlib.crc32(part, crc); x ^= xor_fold(part)
                c[CTL_R1], c[CTL_R2], c[CTL_R3] = s & 0xFF, x, crc; c[CTL_STATUS] = 1
        elif cmd == 8:
            if 0 <= off <= n and 0 <= cnt and off + cnt <= MAX_FILE:
                op = self._prep(off, cnt); v = a3 & 0xFF
                for p in range(off, off + cnt, CHUNK):
                    m = min(CHUNK, off + cnt - p); self.file[p:p + m] = bytes((v,)) * m
                self._commit(op, off, cnt); c[CTL_R1] = cnt; c[CTL_STATUS] = 1
        elif cmd == 9:
            if 0 <= off and 0 <= cnt and off + cnt <= n and 0 <= a3 <= n and a3 + cnt <= MAX_FILE:
                op = self._prep(a3, cnt)
                if a3 <= off or a3 >= off + cnt: starts = range(0, cnt, CHUNK)
                else: starts = range((cnt - 1) // CHUNK * CHUNK, -1, -CHUNK)
                for p in starts:
                    m = min(CHUNK, cnt - p); self.file[a3 + p:a3 + p + m] = self.file[off + p:off + p + m]
                self._commit(op, a3, cnt); c[CTL_R1] = cnt; c[CTL_STATUS] = 1
        elif cmd == 10:
            if 0 <= off and 0 <= cnt and off + cnt <= n:
                op = self._prep(off, cnt); tbl = bytes(i ^ (a3 & 0xFF) for i in range(256))
                for p in range(off, off + cnt, CHUNK):
                    m = min(CHUNK, off + cnt - p); self.file[p:p + m] = bytes(self.file[p:p + m]).translate(tbl)
                self._commit(op, off, cnt); c[CTL_R1] = cnt; c[CTL_STATUS] = 1
        elif cmd == 11:
            if 0 <= off <= MAX_FILE:
                if off < n:
                    removed = bytes(self.file[off:])
                    if len(removed) > UNDO_CAP // 2: self.u_drop()
                    else: self.u_rec(('d', off, removed), 120 + len(removed))
                    del self.file[off:]; self.dirty = 1
                elif off > n: self._prep(n, off - n)
                c[CTL_STATUS] = 1
        else: warn(f"unknown file command {cmd}")
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
            self.push(self.fetch(self.pop())); self.pc+=1
        elif op=="STORE":
            addr,val=self.pop(),self.pop(); self.store(addr,val); self.pc+=1
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
def xor_fold(data):
    """xor of all bytes of `data` using big-int folding (fast for megabyte chunks)"""
    if not data: return 0
    size = 1
    while size < len(data): size <<= 1
    v = int.from_bytes(bytes(data) + bytes(size - len(data)), "big"); w = size * 8
    while w > 8:
        w //= 2; v = (v >> w) ^ (v & ((1 << w) - 1))
    return v
def cmd_opcodes():
    by={}
    for c,n in sorted(CODON_TO_OP.items()): by.setdefault(n,[]).append(c)
    print(f"HolyGATC {VERSION} — {len(OPDEF)} opcodes, {len(CODON_TO_OP)} codons\n")
    print(f"{'Sym':<6} {'Name':<7} {'Byte':<6} {'Codon':<8} Stack")
    print("-"*62)
    for sym,name,byte,imm,stack in OPDEF:
        print(f"{sym:<6} {name:<7} 0x{byte:02X}  {','.join(by.get(name,['-'])):<8} {stack}")
def main(argv=None):
    p=argparse.ArgumentParser(prog="gatcload")
    p.add_argument("command", nargs="?", default="help", choices=["run","compile","disasm","translate","opcodes","help","version"])
    p.add_argument("source", nargs="?")
    p.add_argument("-o","--output"); p.add_argument("--dna", action="store_true")
    p.add_argument("--max-steps", type=int, default=STEP_LIMIT_DEFAULT)
    p.add_argument("--no-banner", action="store_true"); p.add_argument("--gene", action="store_true")
    p.add_argument("--load", help="load a binary file into VM memory at address 0")
    p.add_argument("--file", help="open FILE (up to 1 GiB) in the file-backed address space used by hexeditor.aa")
    p.add_argument("--net", action="store_true", help="allow the program to use the network (ping, dig, whois, TCP/TLS) - needed by netkit.aa")
    p.add_argument("--save", help="where the editor's w command writes (default: the --file path, in place)")
    p.add_argument("--fs", action="store_true", help="allow the program to use the file system (ls, cp, mv, rm, ...) - needed by sak.aa")
    p.add_argument("--fs-root", metavar="DIR", help="with --fs: start in DIR and refuse paths outside it (best-effort jail)")
    args=p.parse_intermixed_args(argv)
    if args.command=="help" and not args.source:
        print(BANNER); p.print_help()
        print("\n  python3 GATCloader.py run program.aa --no-banner")
        print("  python3 GATCloader.py run netkit.aa --no-banner --net --max-steps 500000000")
        print("  python3 GATCloader.py run sak.aa --no-banner --fs --max-steps 2000000000")
        print("  python3 GATCload.py run hex.aa --no-banner --load /etc/passwd")
        print("  python3 GATCload.py compile program.aa --dna -o program.dna")
        return 0
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
        vm = VM(program, genome=genome)
        if args.load:
            vm.load_file(args.load)
        vm.savepath = args.save
        if args.file:
            vm.load_edit_file(args.file)
        if args.net:
            vm.enable_net()
        if args.fs or args.fs_root:
            vm.enable_fs(args.fs_root)
        try:
            return vm.run(max_steps=args.max_steps)
        except BrokenPipeError:            # e.g. piped into head: stop quietly
            try: sys.stdout = open(os.devnull, "w")
            except OSError: pass
            return 0
        finally:
            if vm.net is not None: vm.net.close()
            if vm.fs is not None: vm.fs.close_all()
    except GATCError as e:
        warn(str(e)); return 1
if __name__=="__main__":
    sys.exit(main())
