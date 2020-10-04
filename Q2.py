import library as poly
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
#let our initial guess be 10
g=10
coeff=[1,-3,-7,27,-18]#coefficients are stored in a list
A=poly.roots(g,coeff)#finding the roots
n=len(coeff)
print("The roots of the equation ",end='')
for i in range(n):
    print("("+str(coeff[i])+")"+"x"+str(n-i-1).translate(SUP),end='')
    if(i!=n-1): print("+",end='')
print(" are :")
for i in range(n-1):
    print("x"+str(i+1).translate(SUB)+" = ",A[i])
print("Which is in accordance with the roots of the equation (3,2,1,-3) with an error margin of 1/10"+str(6).translate(SUP))
#Output
'''
The roots of the equation (1)x⁴+(-3)x³+(-7)x²+(27)x¹+(-18)x⁰ are :
x₁ =  3.0000000000000018
x₂ =  1.9999999999999958
x₃ =  1.0000000000000027
x₄ =  -3.0
Which is in accordance with the roots of the equation (3,2,1,-3) with an error margin of 1/10⁶
'''