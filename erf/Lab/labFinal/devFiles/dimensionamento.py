import sympy as sp

# 1. Definir Símbolos
Vcc, Vbe, Vce, Ic, Ib, Ie, beta = sp.symbols('Vcc Vbe Vce Ic Ib Ie beta')
Rb, Vth, R1, R2, Re, Rc = sp.symbols('Rb Vth R1 R2 Re Rc')
Ve, Vc = sp.symbols('Ve Vc')

# 2. Valores que tu escolhes (Exemplo baseado no teu print do LTSpice)
valores_fixos = {
    Vbe: 0.9,
    beta: 72.534,
    Ic: 4e-3,
    Rc: 1350,
    Re: 200,
    Vcc: 10,
    Vce: 5

}

# 3. Equações de Análise (Caminho Inverso)
equacoes = [
    sp.Eq(Ib, Ic/beta),
    sp.Eq(Ie, Ic+Ib),
    sp.Eq(Rb, Re*15),
    sp.Eq(Ve, Ie*Re),
    sp.Eq(Vc, Vcc-Vce-Ve),
    sp.Eq(Vth, Rb * Ib + Vbe + Ve),
    sp.Eq(Vth, Vcc * (R2 / (R1 + R2))),
    sp.Eq(Rb, (R1 * R2) / (R1 + R2))
]

# 4. Resolver para as variáveis de operação (Ic, Vce)
res = sp.solve([eq.subs(valores_fixos) for eq in equacoes], [Ic, Ie, Ib, Rb, Ve, Vc, Vth, R1, R2], dict=True)[0]

print(f"--- Ponto de Operação Resultante ---")
print(f"R1: {float(res[R1]/1000):.1f}k Ohms")
print(f"R2: {float(res[R2]/1000):.1f}k Ohms")