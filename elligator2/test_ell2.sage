
load("ell2generic.sage")
load("ell2curve255.sage")

Hex = lambda x: map(hex,map(int,x))

class SuiteTest(SageObject):
    def _test_curve255(self,tester):
        """
        Generic vs curve25519
        """
        # Curve25519
        p = 2**255 - 19
        A = 486662
        B = 1
        ell2Gen = Ell2Generic(A,B,p)
        ell2cu255 = Ell2Curve255()
        for i in range(100):
            e = randint(0,p)
            P0 = ell2Gen.maptocurve(e)
            P1 = ell2cu255.maptocurve(e)
            tester.assertEqual(P0 , P1)

    # def _test_curve448(self,tester):
    #     """
    #     Generic vs curve448
    #     """
    #     # Curve448
    #     p = 2^448 - 2^224 - 1
    #     A = 156326
    #     B = 1
    #     ell2Gen = Ell2Generic(A,B,p)
    #     ell2cu448 = Ell2Curve448()
    #     for i in range(100):
    #         e = randint(0,p)
    #         P0 = ell2Gen.maptocurve(e)
    #         P1 = ell2cu448.maptocurve(e)
    #         tester.assertEqual(P0 , P1)

TestSuite(SuiteTest()).run(
    verbose = True,
    skip    = ["_test_pickling",
               "_test_new",
               "_test_category",
               "_test_not_implemented_methods"]
)
