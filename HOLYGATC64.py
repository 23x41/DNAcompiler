#!/usr/bin/env python3
"""HOLYGATC64 — stream a real FASTA without loading it.

GRCh38 is gigabytes. The old compiler called strip_to_bases on the whole
file, then compiled every codon, then ran. This one never does that.
It reads bases until the window is full, executes each codon, and moves on.

  python3 HOLYGATC64.py run human.fna --bases 3000
  python3 HOLYGATC64.py run human.fna --skip 100000 --bases 10000
  python3 HOLYGATC64.py run program.aa

Genome rules: HALT and QUIT do not stop. JMP CALL JZ JNZ RET walk forward.
IN and GETC push -1 and do not touch the keyboard. Faults print once.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

try:
    import GATC64 as G
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import GATC64 as G

FASTA = (".fna", ".fa", ".fasta", ".dna", ".gatc")

def warn(msg):
    G.warn(msg)

def peek_packed(path):
    lines = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for _ in range(30):
            line = f.readline()
            if not line:
                break
            lines.append(line)
    blob = "".join(lines)
    return G.looks_packed(blob), blob

def stream_codons(path, skip, limit):
    """Yield codons. Stop after `limit` bases of sequence (0 means no cap)."""
    seen = 0
    buf = []
    forced = False
    with open(path, encoding="utf-8", errors="replace") as f:
        first = True
        for line in f:
            if line.startswith(">"):
                continue
            for ch in line:
                cu = ch.upper()
                if cu == "U":
                    cu = "T"
                if cu not in "ACGT":
                    continue
                if first:
                    first = False
                    if cu != "A":
                        # may not be ATG; decide after 3 bases
                        pass
                if seen < skip:
                    seen += 1
                    continue
                buf.append(cu)
                seen += 1
                if not forced and len(buf) == 3 and "".join(buf) != "ATG":
                    warn("no start codon; forcing ATG and skipping it")
                    forced = True
                    # do not emit the forced codon; keep the real bases
                while len(buf) >= 3:
                    codon = "".join(buf[:3])
                    del buf[:3]
                    if forced == "skip":
                        pass
                    yield codon
                if limit and (seen - skip) >= limit:
                    return
        if buf and not limit:
            warn("trailing bases " + "".join(buf) + "; ignored")

class StreamVM(G.VM):
    def pop(self):
        if not self.stack:
            warn("stack underflow; using 0")
            return 0
        return self.stack.pop()

def exec_codon(vm, codon):
    name = G.CODON_TO_OP.get(codon, "NOP")
    if name == "LIT":
        return "LIT"
    if name in ("IN", "GETC"):
        warn(name + " in genome; pushing -1")
        vm.push(-1)
        return None
    if name in ("JMP", "JZ", "JNZ", "CALL", "RET"):
        warn(name + " in genome; walking forward")
        if name != "RET":
            vm.pop()
            if name in ("JZ", "JNZ"):
                vm.pop()
        return None
    if name in ("HALT", "QUIT"):
        warn(name + " in genome; continuing")
        if name == "QUIT":
            vm.pop()
        return None
    vm._exec(G.Instr(name))
    # _exec advances pc; streaming mode ignores pc
    return None

def run_fasta(path, skip, bases, progress):
    vm = StreamVM([], genome=True)
    n = 0
    pending_lit = False
    print(f"streaming {path} skip={skip} bases={bases or 'all'}", file=sys.stderr)
    for codon in stream_codons(path, skip, bases):
        if pending_lit:
            vm.push(G.codon_value(codon))
            pending_lit = False
            n += 1
            continue
        kind = exec_codon(vm, codon)
        if kind == "LIT":
            pending_lit = True
        n += 1
        if progress and n % progress == 0:
            print(f"codons {n}  stack {len(vm.stack)}", file=sys.stderr)
    if pending_lit:
        warn("LIT at end of window; pushing 0")
        vm.push(0)
    vm.stdout.flush()
    print(f"done {n} codons  stack {vm.stack[-8:]}", file=sys.stderr)
    return 0

def main(argv=None):
    ap = argparse.ArgumentParser(description="HOLYGATC64 streaming genome runner")
    ap.add_argument("command", nargs="?", default="run")
    ap.add_argument("source", nargs="?")
    ap.add_argument("--bases", type=int, default=0, help="stop after this many bases (0 = whole file)")
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--progress", type=int, default=1000)
    ap.add_argument("--no-banner", action="store_true")
    args = ap.parse_args(argv)
    if args.command in ("help", "-h"):
        ap.print_help()
        return 0
    if not args.source:
        warn("source file required")
        return 2
    lower = args.source.lower()
    if not args.no_banner:
        print("HOLYGATC64  stream runner  halt does not stop", file=sys.stderr)
    if lower.endswith(FASTA):
        packed, _ = peek_packed(args.source)
        if packed and args.skip == 0 and args.bases == 0:
            text = open(args.source, encoding="utf-8", errors="replace").read()
            program, genome = G.compile_source(text, args.source)
            return G.VM(program, genome=genome).run()
        # A multi-gigabyte default of "all" will run for hours. Say so, then stream.
        if args.bases == 0:
            print("no --bases given; streaming the whole file, progress every %d codons" % args.progress, file=sys.stderr)
        return run_fasta(args.source, args.skip, args.bases, args.progress)
    text = open(args.source, encoding="utf-8", errors="replace").read()
    program, genome = G.compile_source(text, args.source)
    return G.VM(program, genome=genome).run()

if __name__ == "__main__":
    sys.exit(main())
