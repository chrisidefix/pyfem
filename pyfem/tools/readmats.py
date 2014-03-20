# -*- coding: utf-8 -*-

# Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

import json

p2pMap = {
    u'ρS'  : 'rhoS',
    u'ν'   : 'nu',
    u'φ'   : 'phi',
    u'nf'  : 'nf',
    u'βl'  : 'betL',
    u'βg'  : 'betG',
    u'ρ'   : 'rho',
    u'ρL'  : 'rhoL',
    u'ρG'  : 'rhoG',
    u'RΘg' : 'Rthg',
    u'λd'  : 'lamd',
    u'λw'  : 'lamw',
    u'βd'  : 'betd',
    u'βw'  : 'betw',
    u'β1'  : 'bet1',
    u'β2'  : 'bet2',
    u'α'   : 'alp',
    u'μ'   : 'mu',
    u'τy0' : 'C',
    u'k1'  : 'kn',
    u'ks'  : 'ks',
}
def p2pFunc(p):
    if p in p2pMap: return p2pMap[p]
    return p

#p2pMap.get(p,p)

def readGofemMat(fnkey):
    # read file
    if fnkey[-4:] == '.mat': fnkey = fnkey[:-4]
    fil  = open('%s.mat' % fnkey, 'r')
    Mori = json.load(fil)
    fil.close()
    # parse materials
    M = {}
    for m in Mori:
        M[m['name']] = {}
        for i, p in enumerate(m['prms']):
            M[m['name']][p2pFunc(p)] = m['vals'][i]
    return M

if __name__ == "__main__":
    mats = readGofemMat('cspaper')
    for name, prms in mats.items():
        print name
        print "{"
        for p, v in prms.items():
            print "  %6s: %g," % (p, v)
        print "}"
