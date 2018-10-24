"""
@author Armando Faz

Implementation of Elligator2 mapping by Bernstein, Hamburg, Krasnova, Lange.


"""

class Ell2Generic:
    """
    A textbook implementation of elligator2.
    """
    def __init__(self,A,B,p):
        """
        Instantiates elligator2 mapping function.
        Maps elements in F=GF(p) to points in the elliptic
        curve E/F: By^2=x^3+Ax^2+x
        """
        self.p = p
        self.F = GF(p)
        self.A = self.F(A)
        self.B = self.F(B)
        self.E = EllipticCurve(self.F, [0, self.A, 0, self.B, 0])

    def maptocurve(self,r):
        """
        Map r to a point in E
        Assumptions:
           The returned point can be in any subgroup of E.
           u = 2 if p=5 mod 8
           u =-1 if p=3 mod 4
        """
        if 3 == self.p%4:
            u = self.F(-1)
        if 5 == self.p%8:
            u = self.F(2)
        assert(u.is_square()==false)

        v = -self.A/(1+u*r**2)
        e = legendre_symbol(v**3+self.A*v**2+self.B*v,self.p)
        x = e*v-(1-e)*self.A/self.F(2)
        y = -e*self._sqrt(x**3+self.A*x**2+self.B*x)
        return self.E([x,y])

    def _abs(self,a):
        """
        Returns the positive element of a.
        Assumptions:
            Positive elements are those in the range [  0  , (p-1)/2 ].
            Negative elements are those in the range [(p-1)/2+1, p-1 ].
        """
        lim = (self.p-1)//2
        if ZZ(a) > ZZ(lim):
            a = -a
        return a

    def _sqrt(self,a):
        """
        Returns the square-root of a.
        Assumptions: Since there are two valid roots of a, one must
        cannoncally choose one of them.
        if p=3 mod 4
           sqrt(u) returns the principal square root of u.
        if p=5 mod 8
           sqrt(u) returns the absolute value of the root of u.
        """
        if 5 == self.p%8 :
            return self._sqrt5mod8(a)
        if 3 == self.p%4 :
            return self._sqrt3mod4(a)
        return None

    def _sqrt3mod4(self,a):
        """
        Returns the principal square-root of a.
        Assumption:
            p == 3 mod 4.
        """
        exp = ZZ( (self.p+1)//4 )
        b = a**exp
        return b

    def _sqrt5mod8(self,a):
        """
        Returns the absolute value of the square-root of a.
        Assumption:
            p == 5 mod 8.
        """
        exp = ZZ( (self.p+3)//8 )
        b = a**exp
        if b**2 != a:
            b = b*sqrt(self.F(-1))
        return self._abs(b)
