import sys, getopt


def main(args):
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv[1:])


main(sys.argv[1:])
