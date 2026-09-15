# HolyGATC opcodes — 36

MADE THE COMPILER WITH SUPERGROK

Stack effects use Forth notation: `(before -- after)`. TOS = top of stack.

| Sym | Name  | Byte | Imm | Stack              | Description                    |
|-----|-------|------|-----|--------------------|--------------------------------|
| A   | ADD   | 0x01 |     | (a b -- a+b)       | integer add                    |
| B   | JMP   | 0x02 |     | (addr --)          | jump to address                |
| C   | CALL  | 0x03 |     | (addr --)          | call subroutine                |
| D   | DUP   | 0x04 |     | (a -- a a)         | duplicate TOS                  |
| E   | EQ    | 0x05 |     | (a b -- flag)      | 1 if a == b                    |
| F   | FETCH | 0x06 |     | (addr -- val)      | load memory[addr]              |
| G   | GT    | 0x07 |     | (a b -- flag)      | 1 if a > b                     |
| H   | HALT  | 0x08 |     | (--)               | stop                           |
| I   | IN    | 0x09 |     | (-- n)             | read integer                   |
| J   | JZ    | 0x0A |     | (flag addr --)     | jump if flag == 0              |
| K   | LIT   | 0x0B | yes | (-- n)             | push immediate                 |
| L   | LT    | 0x0C |     | (a b -- flag)      | 1 if a < b                     |
| M   | MUL   | 0x0D |     | (a b -- a*b)       | multiply                       |
| N   | NOT   | 0x0E |     | (a -- flag)        | 1 if a == 0                    |
| O   | OR    | 0x0F |     | (a b -- a\|b)      | bitwise OR                     |
| P   | AND   | 0x10 |     | (a b -- a&b)       | bitwise AND                    |
| Q   | QUIT  | 0x11 |     | (code --)          | exit with status               |
| R   | RET   | 0x12 |     | (--)               | return from CALL               |
| S   | SUB   | 0x13 |     | (a b -- a-b)       | subtract                       |
| T   | STORE | 0x14 |     | (val addr --)      | memory[addr] = val             |
| U   | OVER  | 0x15 |     | (a b -- a b a)     | copy NOS                       |
| V   | DIV   | 0x16 |     | (a b -- a/b)       | integer divide                 |
| W   | PUTC  | 0x17 |     | (n --)             | write character                |
| X   | XOR   | 0x18 |     | (a b -- a^b)       | bitwise XOR                    |
| Y   | PUTN  | 0x19 |     | (n --)             | write decimal                  |
| Z   | SWAP  | 0x1A |     | (a b -- b a)       | swap                           |
| _   | NOP   | 0x1B |     | (--)               | no operation                   |
| %   | MOD   | 0x1C |     | (a b -- a%b)       | modulo                         |
| ~   | NEG   | 0x1D |     | (a -- -a)          | negate                         |
| !   | DROP  | 0x1E |     | (a --)             | discard TOS                    |
| ?   | GETC  | 0x1F |     | (-- n)             | read character                 |
| <   | SHL   | 0x20 |     | (a b -- a<<b)      | shift left                     |
| >   | SHR   | 0x21 |     | (a b -- a>>b)      | shift right                    |
| $   | RAND  | 0x22 |     | (n -- r)           | random in [0, n)               |
| +   | INC   | 0x23 |     | (a -- a+1)         | increment                      |
| -   | DEC   | 0x24 |     | (a -- a-1)         | decrement                      |

Word aliases also accepted: `ADD`, `PUSH`, `LIT`, `HALT`, `PRINT`, `EMIT`, `LOAD`, `STORE`, …
