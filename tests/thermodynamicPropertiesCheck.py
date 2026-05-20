import math

def JanafCoeffs(TTR, specie):
    if (specie == "N2"):
        lowCoeffs = [  2.210371497e+04, -3.818461820e+02, 6.082738360e+00,
                      -8.530914410e-03,  1.384646189e-05,-9.625793620e-09,
                       2.519705809e-12,  7.108460860e+02,-1.076003744e+01]
        midCoeffs = [  5.877124060e+05, -2.239249073e+03, 6.066949220e+00,
                      -6.139685500e-04,  1.491806679e-07,-1.923105485e-11,
                       1.061954386e-15,  1.283210415e+04,-1.586640027e+01]
        highCoeffs = [ 8.310139160e+08, -6.420733540e+05, 2.020264635e+02,
                      -3.065092046e-02,  2.486903333e-06,-9.705954110e-11,
                       1.437538881e-15,  4.938707040e+06,-1.672099740e+03]
    elif (specie == "O2"):
        lowCoeffs = [  -3.425563420e+04, 4.847000970e+02, 1.119010961e+00,
                        4.293889240e-03,-6.836300520e-07,-2.023372700e-09,
                        1.039040018e-12,-3.391454870e+03, 1.849699470e+01]
        midCoeffs = [  -1.037939022e+06, 2.344830282e+03, 1.819732036e+00,
                        1.267847582e-03,-2.188067988e-07, 2.053719572e-11,
                       -8.193467050e-16,-1.689010929e+04, 1.738716506e+01]
        highCoeffs = [ 4.975294300e+08, -2.866106874e+05, 6.690352250e+01,
                      -6.169959020e-03,  3.016396027e-07,-7.421416600e-12,
                       7.278175770e-17,  2.293554027e+06,-5.530621610e+02]
    elif (specie == "N"):
        lowCoeffs = [  0.000000000e+00,  0.000000000e+00, 2.500000000e+00,
                       0.000000000e+00,  0.000000000e+00, 0.000000000e+00,
                       0.000000000e+00,  5.610463780e+04, 4.193905036e+00]
        midCoeffs = [   8.876501380e+04,-1.071231500e+02, 2.362188287e+00,
                        2.916720081e-04,-1.729515100e-07, 4.012657880e-11,
                       -2.677227571e-15, 5.697351330e+04, 4.865231506e+00]
        highCoeffs = [ 5.475181050e+08, -3.107574980e+05, 6.916782740e+01,
                      -6.847988130e-03,  3.827572400e-07,-1.098367709e-11,
                       1.277986024e-16,  2.550585618e+06,-5.848769753e+02]
    elif (specie == "O"):
        lowCoeffs = [  -7.953611300e+03, 1.607177787e+02, 1.966226438e+00,
                        1.013670310e-03,-1.110415423e-06, 6.517507500e-10,
                       -1.584779251e-13, 2.840362437e+04, 8.404241820e+00]
        midCoeffs = [   2.619020262e+05,-7.298722030e+02, 3.317177270e+00,
                       -4.281334360e-04, 1.036104594e-07,-9.438304330e-12,
                        2.725038297e-16, 3.392428060e+04,-6.679585350e-01]
        highCoeffs = [ 1.779004264e+08, -1.082328257e+05, 2.810778365e+01,
                      -2.975232262e-03,  1.854997534e-07,-5.796231540e-12,
                       7.191720164e-17,  8.890942630e+05,-2.181728151e+02]
    elif (specie == "NO"):
        lowCoeffs = [  -1.143916503e+04,  1.536467592e+02, 3.431468730e+00,
                       -2.668592368e-03,  8.481399120e-06, -7.685111050e-09,
                        2.386797655e-12,  9.098214410e+03,  6.728725490e+00]
        midCoeffs = [   2.239018716e+05, -1.289651623e+03,  5.433936030e+00,
                       -3.656034900e-04,  9.880966450e-08, -1.416076856e-11,
                        9.380184620e-16,  1.750317656e+04, -8.501669090e+00]
        highCoeffs = [ -9.575303540e+08,  5.912434480e+05, -1.384566826e+02,
                        1.694339403e-02, -1.007351096e-06,  2.912584076e-11,
                       -3.295109350e-16, -4.677501240e+06,  1.242081216e+03]
    if (TTR <= 1000):
        return lowCoeffs
    elif (TTR > 1000 and TTR <= 6000):
        return midCoeffs
    else:
        return highCoeffs

def JanafCp(coeffs, TTR):
    a = coeffs
    Cp = a[0]/(TTR*TTR) + a[1]/TTR + a[2] + (a[3] + (a[4] + (a[5] + a[6]*TTR)*TTR)*TTR)*TTR
    return Cp

def CvT(Cp, R_s):
    Cv = Cp - R_s
    return Cv
    
def CvR(R_s):
    return R_s
    
def CvVib(TTR, thetaVib, R_s):
    CvV = R_s*(thetaVib/(2*TTR))/(math.sinh(thetaVib/(2*TTR)))
    return CvV
    
def mu(TTR, specie):
    if (specie == "N2"):
        AB = 2.7321e-02
        BB = 2.4321e-01
        CB = -1.0805e+01
    elif (specie == "O2"):
        AB = 2.6757e-02
        BB = 2.5401e-01
        CB = -1.0711e+01
    elif (specie == "N"):
        AB = 2.4768e-02
        BB = 2.8672e-01
        CB = -1.0976e+01
    elif (specie == "O"):
        AB = 2.7230e-02
        BB = 2.3976e-01
        CB = -1.0383e+01
    elif (specie == "NO"):
        AB = 2.7227e-02
        BB = 2.4463e-01
        CB = -1.0776e+01
    mu = 0.1*math.exp((AB*math.log(TTR) + BB)*math.log(TTR) + CB)
    return mu

def kappaTR(TTR, specie, R_s, Cp):
    kappaT = mu(TTR, specie) * (2.5 * CvT(Cp, R_s) + CvR(R_s))
    return kappaT
    
def kappaVib(TTR, specie):
    nu = 1.2
    kappaV = nu * mu(TTR, specie)
    return kappaV
    
def printValues(specie, R_s, CvT, CvR, CpTR, CvVib, mu, kappaTR, kappaVib):
    if (specie == "N2" or specie == "O2" or specie == "NO"):
        print("***********************************************************************************")
        print(f"For {specie}, calculated thermodynamical values are given below")
        print("Specific gas constant is: ", R_s)
        print("Average Translational Cv value in internal cells should be ", CvT)
        print("Average Rotational Cv value in internal cells should be ", CvR)
        print("Average Translational-Rotational Cp value in internal cells should be ", CpTR)
        print("Average Vibrational Cv value in internal cells should be ", CvVib)
        print("Average mu value in internal cells should be ", mu)
        print("Average Translational-Rotational kappa value in internal cells should be ", kappaTR)
        print("Average Vibrational kappa value in internal cells should be ", kappaVib)
    else:
        print("***********************************************************************************")
        print(f"For {specie}, calculated thermodynamical values are given below")
        print("Specific gas constant is: ", R_s)
        print("Average Translational Cv value in internal cells should be ", CvT)
        print("Average Rotational Cv value in internal cells should be ", CvR)
        print("Average Translational-Rotational Cp value in internal cells should be ", CpTR)
        print("Average mu value in internal cells should be ", mu)
        print("Average Translational-Rotational kappa value in internal cells should be ", kappaTR)
        
print("Sanity Check for neTCLib Implementation")

print(" ")

R = 8.314

specie_N2 = "N2"
specie_O2 = "O2"
specie_N  = "N"
specie_O  = "O"
specie_NO = "NO"
molWeight_N2 = 28.01348e-03
molWeight_O2 = 31.9988e-03
molWeight_N = 14.0067e-03
molWeight_O = 15.9994e-03
molWeight_NO = 30.0061e-03
thetaVib_N2 = 3371
thetaVib_O2 = 2256
thetaVib_NO = 2719
R_N2 = R / molWeight_N2
R_O2 = R / molWeight_O2
R_N = R / molWeight_N
R_O = R / molWeight_O
R_NO = R / molWeight_NO
TTR       = 2000
TVib      = 1000

coeffs_N2    = JanafCoeffs(TTR, specie_N2)
CpTR_N2      = JanafCp(coeffs_N2, TTR) * R_N2
CvT_N2       = CvT(CpTR_N2, R_N2)
CvR_N2       = CvR(R_N2)
CvVib_N2     = CvVib(TVib, thetaVib_N2, R_N2)
mu_N2        = mu(TTR, specie_N2)
kappaTR_N2   = kappaTR(TTR, specie_N2, R_N2, CpTR_N2)
kappaVib_N2  = kappaVib(TTR, specie_N2)

coeffs_O2    = JanafCoeffs(TTR, specie_O2)
CpTR_O2      = JanafCp(coeffs_O2, TTR) * R_O2
CvT_O2       = CvT(CpTR_O2, R_O2)
CvR_O2       = CvR(R_O2)
CvVib_O2     = CvVib(TVib, thetaVib_O2, R_O2)
mu_O2        = mu(TTR, specie_O2)
kappaTR_O2   = kappaTR(TTR, specie_O2, R_O2, CpTR_O2)
kappaVib_O2  = kappaVib(TTR, specie_O2)

coeffs_N    = JanafCoeffs(TTR, specie_N)
CpTR_N      = JanafCp(coeffs_N, TTR) * R_N
CvT_N       = CvT(CpTR_N, R_N)
CvR_N       = CvR(R_N)
mu_N        = mu(TTR, specie_N)
kappaTR_N   = kappaTR(TTR, specie_N, R_N, CpTR_N)

coeffs_O    = JanafCoeffs(TTR, specie_O)
CpTR_O      = JanafCp(coeffs_O, TTR) * R_O
CvT_O       = CvT(CpTR_O, R_O)
CvR_O       = CvR(R_O)
mu_O        = mu(TTR, specie_O)
kappaTR_O   = kappaTR(TTR, specie_O, R_O, CpTR_O)

coeffs_NO    = JanafCoeffs(TTR, specie_NO)
CpTR_NO      = JanafCp(coeffs_NO, TTR) * R_NO
CvT_NO       = CvT(CpTR_NO, R_NO)
CvR_NO       = CvR(R_NO)
CvVib_NO     = CvVib(TVib, thetaVib_NO, R_NO)
mu_NO        = mu(TTR, specie_NO)
kappaTR_NO   = kappaTR(TTR, specie_NO, R_NO, CpTR_NO)
kappaVib_NO  = kappaVib(TTR, specie_NO)

printValues(specie_N2, R_N2, CvT_N2, CvR_N2, CpTR_N2, CvVib_N2, mu_N2, kappaTR_N2, kappaVib_N2)
printValues(specie_O2, R_O2, CvT_O2, CvR_O2, CpTR_O2, CvVib_O2, mu_O2, kappaTR_O2, kappaVib_O2)
printValues(specie_N, R_N, CvT_N, CvR_N, CpTR_N, 0, mu_N, kappaTR_N, 0)
printValues(specie_O, R_O, CvT_O, CvR_O, CpTR_O, 0, mu_O, kappaTR_O, 0)
printValues(specie_NO, R_NO, CvT_NO, CvR_NO, CpTR_NO, CvVib_NO, mu_NO, kappaTR_NO, kappaVib_NO)

print("Check is finished!")
