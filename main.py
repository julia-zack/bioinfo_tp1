"""Entry point for the TP. Runs one exercise at a time, by name:

    python3 main.py 4a
    python3 main.py 1b-ii
    python3 main.py          # lists the available exercises

The names match the numbering in the TP PDF, so the command you type, the
section in the report and the file on disk all agree.
"""

import sys

from exercises import ex1, ex2, ex3
from exercises.ex4 import a as ex4a, b as ex4b, c as ex4c, d as ex4d, e as ex4e

EXERCISES = {
    "1a-i":   ex1.a_i,
    "1a-ii":  ex1.a_ii,
    "1b-i":   ex1.b_i,
    "1b-ii":  ex1.b_ii,
    "1b-iii": ex1.b_iii,

    "2a": ex2.a,
    "2b": ex2.b,

    "3b": ex3.b,
    "3c": ex3.c,

    "4a": ex4a.run,
    "4b": ex4b.run,
    "4c": ex4c.run,
    "4d": ex4d.run,
    "4e": ex4e.run,
}


def usage():
    print("usage: python3 main.py <exercise>")
    print("available exercises:")
    for name in EXERCISES:
        print(f"  {name}")


def main():
    if len(sys.argv) != 2:
        usage()
        return
    name = sys.argv[1]
    if name not in EXERCISES:
        print(f"'{name}' is not a known exercise.\n")
        usage()
        return
    EXERCISES[name]()


if __name__ == "__main__":
    main()
