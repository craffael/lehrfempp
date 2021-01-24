#!/usr/bin/python

from __future__ import print_function


def main():
    from argparse import ArgumentParser
    from platform import system
    from subprocess import Popen, PIPE
    from sys import stderr

    parser = ArgumentParser()

    parser.add_argument('args', nargs='*', help='arguments passed to doxygen')
    parser.add_argument(
        '--exe',
        default='doxygen.exe' if system() == 'Windows' else 'doxygen',
        help='path to doxygen executable, default: %(default)s')

    args = parser.parse_args()

    process = Popen([args.exe] + args.args, stdout=PIPE, stderr=PIPE)
    errors = process.stderr.read().strip()
    if len(errors) != 0:
        print(errors, file=stderr)
        exit(2)


try:
    main()
except KeyboardInterrupt:
    pass