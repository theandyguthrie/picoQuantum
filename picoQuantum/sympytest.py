from sympy import MatrixSymbol, Matrix, symbols

a00, a01, a10, a11 = symbols('a00 a01 a10 a11')
b00, b01, b10, b11 = symbols('a00 a01 a10 a11')
A = Matrix([[a00, a01],[a10, a11]])
B = Matrix([[b00,b01],[b10, b11]])


C = A*B
C = C.subs(a00, 5)


print(C)