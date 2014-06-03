#!/usr/bin/env python

import sys

FILE0 = "../history_base.txt"

STEP = 20


def usage():
    print 'usage: %s master_file file' % sys.argv[0]


def read_vals(file):
    f = open(file, 'r')
    vals = {}
    for line in f:
        field = line.split()
        if field[0] == 'step': continue
        if int(field[0]) == STEP:
            vals['v_max']  = field[2]
            vals['v_div_max']  = field[4]
            vals['r_b']    = field[6]
            vals['deltaP'] = field[7]
            vals['avrP']   = field[8]
            vals['deltaV'] = field[9]
            vals['avrV']   = field[10]
            return vals
    return vals


def check(key, vals0, vals):
    if vals0[key] != vals[key]:
        print '%s: %s != %s' % (key, vals0[key], vals[key])
        return 1
    return 0


if __name__ == '__main__':

    if len(sys.argv) != 3:
        usage()
        sys.exit(1)

    print
    print 'Checking result ...'

    vals0 = read_vals(sys.argv[1])
    vals  = read_vals(sys.argv[2])

    ndiff = 0
    ndiff = ndiff + check('v_max', vals0, vals)
   #ndiff = ndiff + check('v_div_max', vals0, vals)
   #ndiff = ndiff + check('r_b', vals0, vals)
   #ndiff = ndiff + check('deltaP', vals0, vals)
   #ndiff = ndiff + check('avrP', vals0, vals)
    ndiff = ndiff + check('deltaV', vals0, vals)
    ndiff = ndiff + check('avrV', vals0, vals)

    if ndiff == 0:
        print 'OK.'
        sys.exit(0)
    else:
        print 'NG.'
        sys.exit(1)

