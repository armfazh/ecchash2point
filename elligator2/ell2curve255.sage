"""
@author Armando Faz

Implementation of Elligator2 for Curve25519.

"""

class Ell2Curve255:
    """
    Optimized implementation of Elligator2 for Curve25519
    """

    def __init__(self,argOptimized=True,argProyective=False):
        """
        Optimized (true default)
        """
        self.p = 2**255-19
        self.F = GF(self.p)
        self.A = self.F(486662)
        self.B = self.F(1)
        self.u = self.F(2)
        self.E = EllipticCurve(self.F, [0, self.A, 0, self.B, 0])
        self.opt = argOptimized
        self.proy = argProyective

    def _abs(self,a):
        """
        Returns the positive element of a.
        Def. Positive elements are those in the range [0 , (p-1)/2 ] .
        """
        lim = (self.p-1)//2
        if ZZ(a) > ZZ(lim):
            a = -a
        return a


    def _sqrt5mod8(self,a):
        """
        Returns the square-root of a for primes congruent to 5 mod 8.
        """
        exp = ZZ( (self.p+3)//8 )
        b = a**exp
        if b**2 != a:
            b = b*sqrt(self.F(-1))
        return self._abs(b)

    def maptocurve(self,r):
        """
        Maps r to a point on Curve25519
        """
        if self.proy:
            return self._proyective(r)
        if self.opt == False:
            return self._ell2curve255(r)
        else:
            return self._expell2curve255(r)

    def _ell2curve255(self,r):
        sqrt_of_minusone = self.F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
        u = self.u
        v = -self.A/(1+u*r**2)
        l = v**3+self.A*v**2+self.B*v
        power = l**ZZ(((self.p+3)//8))
        if power**4 == l**2 :
            # print("case e=+1\n")
            x = v
            y = power
            if power**2 != l:
                # print("case e=+1 neg\n")
                y = y*sqrt_of_minusone
            y = -self._abs(y)
        else:
            # print("case e=-1\n")
            x = v*u*r**2
            l = l*u*r**2
            y = power*r*u**ZZ((self.p+3)//8)
            if y**2 != l:
                # print("case e=-1 neg\n")
                y = y*sqrt_of_minusone
            y = self._abs(y)
        return self.E([x,y])

    def _expell2curve255(self,r):
        sqrt_of_minusone = self.F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
        u = self.u
        v0 = -self.A
        v1 = (1+u*r**2)
        #l = v**3+A*v**2+B*v
        l0 = v0**3+self.A*v0**2*v1+self.B*v0*v1**2
        l1 = v1**3
        #power = l**ZZ(((p+3)//8))
        power = l0*l1**3*(l0*l1**7)**ZZ((self.p-5)//8)
        if l1**2*power**4 == l0**2 :
            # print("case e=+1 pos\n")
            # x = v
            x0 = v0
            x1 = v1
            y = power
            if l1*power**2 != l0:
                # print("case e=+1 neg\n")
                y = y*sqrt_of_minusone
            y = -self._abs(y)
        else:
            # print("case e=-1 pos\n")
            # x = v*u*r**2
            x0 = v0*u*r**2
            x1 = v1
            l0 = l0*u*r**2
            y = power*r*u**ZZ((self.p+3)//8)
            if l1*y**2 != l0:
                # print("case e=-1 neg\n")
                y = y*sqrt_of_minusone
            y = self._abs(y)
        return self.E([x0/x1,y])

    def _proyective(self,r):
        sqrtnegone = self.F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
        ccA3 = self.A**3
        negA = -self.A
        pminus5div8 = ZZ((self.p-5)//8)
        twopowpminus3div8 = self.F(2)**ZZ((self.p+3)//8)
                                      #  M | S | A | M_A | E
        ur2 = 2*r**2                  #    | 1 | 1 |     |
        Z = ur2 + 1                   #    |   | 1 |     |
        d = (ccA3-self.A*Z)*Z-ccA3    #  1 |   | 2 |  1  |
        d2 = d**2                     #    | 1 |   |     |
        Z3 = Z**2*Z                   #  1 | 1 |   |     |
        Z6 = Z3**2                    #    | 1 |   |     |
        Z9 = Z6*Z3                    #  1 |   |   |     |
        Zc = Z6**2                    #    | 1 |   |     |
        Y = d*Z9                      #  1 |   |   |     |
        Y = Y*(Y*Zc)**pminus5div8     #  2 |   |   |     | 1
        Y2 = Y**2                     #    | 1 |   |     |
        Y4 = Y2**2                    #    | 1 |   |     |
        bit1 = Z6*Y4 == d2            #  1 |   |   |     |
        ur2 = 1 if bit1 else ur2      #    |   |   |     |
        X = negA*ur2                  #    |   |   |  1  |
        # Ycoordinate
        d = d*ur2                     #  1 |   |   |     |
        Y = Y * (1 if bit1 else r*twopowpminus3div8) # 1M_A
        Y2 = Y**2                     #    | 1 |   |     |
        bit2 = Z3*Y2 != d             #  1 |   |   |     |
        Y = Y * (sqrtnegone if bit2 else 1) # 1M
        Y = self._abs(Y) * (-1 if bit1 else 1) # 1A
        Y = Y*Z                       #  1 |   |   |     |
        return self.E([X,Y,Z])  # Cost: 11M+8S+5A+3M_A+1E
