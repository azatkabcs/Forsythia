"""
Class Montgomery curves
"""
from gfp2 import GFp2element
from gfp2 import p


class MontyCurve:
    A = None  # GFp2element(0, 0, 0, 16)
    B = None  # GFp2element(0, 0, 0, 16)
    C = None  # GFp2element(0, 0, 0, 16)
    Ap24 = None
    Am24 = None
    C24 = None
    ap24 = None

    def __init__(self, A, B, C):
        if isinstance(A, GFp2element):
            self.A = A
        else:
            self.A = GFp2element(A, 0, 16)
        if isinstance(B, GFp2element):
            self.B = B
        else:
            self.B = GFp2element(B, 0, 16)
        if isinstance(C, GFp2element):
            self.C = C
        else:
            self.C = GFp2element(C, 0, 16)
        self.Ap24 = self.A + self.C * 2
        self.Am24 = self.A - self.C * 2
        self.C24 = self.C * 4
        self.ap24 = self.Ap24 / self.C

    def __str__(self):
        return 'y^2 = ' + "x^3 + (" + str(self.A / self.C) + ") * x^2 + x"

    def jinv(self):
        j = self.A * self.A
        t1 = self.C * self.C
        t0 = t1 + t1
        t0 = j - t0
        t0 = t0 - t1
        j = t0 - t1
        t1 = t1 * t1
        j = j * t1
        t0 = t0 + t0
        t0 = t0 + t0
        t1 = t0 * t0
        t0 = t1 * t0
        t0 = t0 + t0
        t0 = t0 + t0
        j = j.modinv()
        j = t0 * j
        return j


def get_a(P, Q, D):
    pass

class MontyPoint:
    X = GFp2element(0, 0, 0)
    #    Y = GFp2element(0, 0, 0)
    Z = GFp2element(0, 0, 0)
    parent = None

    def __init__(self, X, Z, parent):
        self.parent = parent
        self.X = X
        #        self.Y = Y
        if not Z is None:
            self.Z = Z
        else:
            self.Z = GFp2element(1, 0, 16)

    def __str__(self):
        return ('(' + str(self.X) + ':' + str(self.Z) + ')')

    def __add__(self, other):
        pass

    def mul2(self):
        t0 = self.X - self.Z
        t1 = self.X + self.Z
        t0 = t0 * t0
        t1 = t1 * t1
        z = t0 * self.parent.C24
        x = z * t1
        t1 = t1 - t0
        t0 = t1 * self.parent.Ap24
        z = z + t0
        z = z * t1
        return MontyPoint(x, z, self.parent)

    def mul3(self):
        t0 = self.X - self.Z
        t2 = t0 * t0
        t1 = self.X + self.Z
        t3 = t1 * t1
        t4 = t1 + t0
        t0 = t1 - t0
        t1 = t4 * t4
        t1 = t1 - t3
        t1 = t1 - t2
        t5 = t3 * self.parent.Ap24
        t3 = t5 * t3
        t6 = t2 * self.parent.Am24
        t2 = t2 * t6
        t3 = t2 - t3
        t2 = t5 - t6
        t1 = t2 * t1
        t2 = t3 + t1
        t2 = t2 * t2
        x = t2 * t4
        t1 = t3 - t1
        t1 = t1 * t1
        z = t1 * t0
        return MontyPoint(x, z, self.parent)

    def mul2e(self, e):
        res = MontyPoint(self.X, self.Z, self.parent)
        for i in range(1, e):
            res = res.mul2()
        return res

    def mul3e(self, e):
        res = MontyPoint(self.X, self.Z, self.parent)
        for i in range(1, e):
            res = res.mul3()
        return res

    def xdbladd(self, P, Q, D):
        assert(P.parent == Q.parent == D.parent)
        t0 = P.X + P.Z
        t1 = P.X - P.Z
        x2p = t0 * t0
        t2 = Q.X - Q.Z
        xpq = Q.X + Q.Z
        t0 = t0 * t2
        z2p = t1 * t1
        t1 = t1 * xpq
        t2 = x2p - z2p
        x2p = x2p * z2p
        xpq = self.parent.ap24 * t2
        zpq = t0 - t1
        z2p = xpq + z2p
        xpq = t0 + t1
        z2p = z2p * t2
        zpq = zpq * zpq
        xpq = xpq * xpq
        zpq = D.X * zpq
        xpq = D.Z * xpq
        return [MontyPoint(x2p, z2p, P.parent), MontyPoint(xpq, zpq, P.parent)]

    def ladder(self, other, e):
        pass

    def diffadd(self, other, diff):
        pass

    def scalar(self, k):
        pass

    def __mul__(self, other):
        return self.scalar(other)
