#!/usr/bin/env python3
"""HolyGATC debugger for .aa and .dna programs.
Put this next to GATC64.py.

  python3 gatc_debug.py human.dna
  python3 gatc_debug.py program.aa
  python3 gatc_debug.py program.dna --trace 30

(dbg) commands:
  s / enter   step one instruction
  n           step over a CALL
  c           run until breakpoint, halt, or end
  b 12        breakpoint at instruction 12
  b           list breakpoints
  d 12        delete breakpoint
  l           list around pc
  p           stacks
  m 0         16 memory cells from 0
  r           reset
  q           quit
Genome runs never read the keyboard. GETC and IN push -1.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
try:
    import GATC64 as G
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import GATC64 as G

class DebugVM(G.VM):
    def __init__(self, program, genome=False):
        super().__init__(program, genome=genome)
        self.breaks = set()
        self.last_fault = ""
    def push(self, n):
        if len(self.stack) >= G.STACK_LIMIT:
            self.last_fault = "stack overflow"
            return
        self.stack.append(int(n))
    def pop(self):
        if not self.stack:
            self.last_fault = "stack underflow; using 0"
            return 0
        return self.stack.pop()
    def step(self):
        n = len(self.program)
        if self.halted or self.pc < 0 or self.pc >= n:
            self.halted = True
            return False
        ins = self.program[self.pc]
        if self.genome and ins.name in ("IN", "GETC"):
            self.push(-1)
            self.pc += 1
            self.last_fault = ins.name + " skipped in genome; pushed -1"
        else:
            self._exec(ins)
        self.steps += 1
        return not self.halted

def load(path, gene=False):
    text = open(path, encoding="utf-8", errors="replace").read()
    if gene:
        return G.compile_gene(text), True
    if path.lower().endswith(".gbc"):
        return G.decode_bytecode(open(path, "rb").read()), False
    return G.compile_source(text, path)

def show(vm, window=7):
    ins = vm.program[vm.pc] if 0 <= vm.pc < len(vm.program) else None
    here = (ins.name + (f" {ins.imm}" if ins.imm is not None else "")) if ins else "END"
    print(f"pc {vm.pc:04d}  {here}   steps {vm.steps}   {'GENOME' if vm.genome else 'PACKED'}")
    print("stack ", vm.stack[-12:])
    print("return", vm.calls[-8:])
    if vm.last_fault:
        print("fault ", vm.last_fault)
    lo, hi = max(0, vm.pc - 2), min(len(vm.program), vm.pc + window)
    for i in range(lo, hi):
        mark = ">>" if i == vm.pc else "  "
        br = "*" if i in vm.breaks else " "
        item = vm.program[i]
        extra = f" {item.imm}" if item.imm is not None else ""
        print(f"{mark}{br} {i:04d}: {item.name}{extra}")

def trace(path, n, gene=False):
    program, genome = load(path, gene)
    vm = DebugVM(program, genome)
    print(f"{path}: {len(program)} instructions, genome={genome}")
    for _ in range(n):
        if vm.halted or not (0 <= vm.pc < len(program)):
            break
        show(vm, window=1)
        vm.step()
        print()
    print(f"stopped at pc {vm.pc} after {vm.steps} steps, halt={vm.halted}")

def repl(path, gene=False):
    program, genome = load(path, gene)
    vm = DebugVM(program, genome)
    print(f"{path}: {len(program)} instructions")
    show(vm)
    while True:
        try:
            line = input("(dbg) ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if line in ("", "s", "step"):
            if not vm.step():
                print("halted")
            show(vm)
        elif line in ("n", "next"):
            if 0 <= vm.pc < len(vm.program) and vm.program[vm.pc].name == "CALL" and not vm.genome:
                depth = len(vm.calls)
                vm.step()
                guard = 0
                while not vm.halted and len(vm.calls) > depth and vm.pc not in vm.breaks and guard < 100000:
                    vm.step(); guard += 1
            else:
                vm.step()
            show(vm)
        elif line in ("c", "cont", "continue"):
            guard = 0
            vm.step()
            while not vm.halted and vm.pc not in vm.breaks and guard < G.STEP_LIMIT_DEFAULT:
                vm.step(); guard += 1
            show(vm)
        elif line.startswith("b"):
            rest = line[1:].strip()
            if not rest:
                print("breaks", sorted(vm.breaks) or "none")
            else:
                vm.breaks.add(int(rest, 0)); print("break", rest)
        elif line.startswith("d") and line[1:].strip():
            vm.breaks.discard(int(line[1:].strip(), 0))
        elif line in ("l", "list"):
            show(vm, window=12)
        elif line in ("p", "stack"):
            print(vm.stack); print("return", vm.calls)
        elif line.startswith("m"):
            addr = int(line[1:].strip() or "0", 0)
            print(addr, vm.mem[addr:addr+16])
        elif line in ("r", "reset"):
            vm = DebugVM(program, genome); show(vm)
        elif line in ("q", "quit"):
            return 0
        else:
            print("s step | n next | c continue | b addr | d addr | l list | p stack | m addr | r reset | q quit")

def main(argv=None):
    ap = argparse.ArgumentParser(description="HolyGATC .dna debugger")
    ap.add_argument("source")
    ap.add_argument("--trace", type=int, default=0)
    ap.add_argument("--gene", action="store_true")
    args = ap.parse_args(argv)
    if args.trace:
        trace(args.source, args.trace, args.gene); return 0
    return repl(args.source, args.gene)

if __name__ == "__main__":
    sys.exit(main())
