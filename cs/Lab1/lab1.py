import math
import random
import matplotlib.pyplot as plt
import numpy as np

Nbits = 12
Vref = 1.0
Vlsb = (2 * Vref) / (2**Nbits)
codigos = []

SC = 0e-2
Soffset = 0e-3
Cu = 1.0
Cpb = 0
Cpl = 0
Cpm = 0

multiplicadores_split_MSB = [16, 8, 4, 2, 1]
multiplicadores_split_C6 = [1]
multiplicadores_split_LSB = [16, 8, 4, 2, 1]

Nmsb = len(multiplicadores_split_MSB)
C6 = len(multiplicadores_split_C6)
Nlsb = len(multiplicadores_split_LSB)

def randn():
    return random.gauss(0, 1)

valores_Cia = [m * (Cu + SC * randn())/2 for m in multiplicadores_split_MSB]
valores_Cib = [m * (Cu + SC * randn())/2 for m in multiplicadores_split_MSB]
valores_CiMSB = [valores_Cia[i] + valores_Cib[i] for i in range(len(valores_Cia))]
valores_CiLSB = [m * (Cu + SC * randn()) for m in multiplicadores_split_LSB]

C6 = 1.0 * Cu
Cb = 2*Nlsb / (2*Nlsb - 1) * Cu + Cpb
Cdl = (2*Nlsb - 1)(Cb-Cu) - Cpl

CtotalMSB = sum(valores_Cia) + sum(valores_Cib) + Cpm
CtotalLSB = sum(valores_CiLSB) + Cdl + Cpl

C_ramo_LSB = (Cb * (CtotalLSB)) / (Cb + CtotalLSB)
Ctot = CtotalMSB + C6 + C_ramo_LSB


# Simulação do processo de conversão
amostras = 8000
Vin_inicial = np.linspace(-Vref, Vref, amostras)



"""def sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref):
    bit = []
    if (Vxp - Vxn) < 0:
        bit.append(0)
    else:
        bit.append(1)

    for j in range(Nmsb):
        if bit[j-1] == 1:
            Vxp -= (valores_CiMSB[j-1] / Ctot) * Vref/2
            Vxn += (valores_CiMSB[j-1] / Ctot) * Vref/2
            if Vxp-Vxn > 0:
                bit.append(1)
            else:
                bit.append(0)
            
        elif bit[j-1] == 0:
            Vxp += (valores_CiMSB[j-1] / Ctot) * Vref/2
            Vxn -= (valores_CiMSB[j-1] / Ctot) * Vref/2

            if Vxp-Vxn > 0:
                bit.append(1)
            else:
                bit.append(0)

    # C6
    if bit[5] == 1:
        Vxp -= (C6 / Ctot) * Vref/2
        Vxn += (C6 / Ctot) * Vref/2
        
        if Vxp-Vxn > 0:
            bit.append(1)
        else:
            bit.append(0)

    elif bit[5] == 0:
        Vxp += (C6 / Ctot) * Vref/2
        Vxn -= (C6 / Ctot) * Vref/2

        if Vxp-Vxn > 0:
            bit.append(1)
        else:
            bit.append(0)

    atenuacao = Cb / (CtotalLSB + Cpl + Cb)
        
    for j in range(Nmsb+2,Nmsb+2+Nlsb):
        if bit[j-1] == 1:
            Vxp -= (valores_CiLSB[j-Nmsb-2] / Ctot) * Vref/2 * atenuacao
            Vxn += (valores_CiLSB[j-Nmsb-2] / Ctot) * Vref/2 * atenuacao

            if Vxp-Vxn > 0:
                bit.append(1)
            else:
                bit.append(0)

        elif bit[j-1] == 0:
            Vxp += (valores_CiLSB[j-Nmsb-2] / Ctot) * Vref/2 * atenuacao
            Vxn -= (valores_CiLSB[j-Nmsb-2] / Ctot) * Vref/2 * atenuacao

            if Vxp-Vxn > 0:
                bit.append(1)
            else:
                bit.append(0)
    return bit"""

def sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref):
    bit = []

    for j in range(Nmsb):
        # Usamos o ÚLTIMO bit decidido para ajustar a tensão
        v_step = (valores_CiMSB[j] / Ctot) * (Vref / 2)
        if bit[-1] == 1:
            Vxp -= v_step
            Vxn += v_step
        else:
            Vxp += v_step
            Vxn -= v_step
        
        # Agora decidimos o PRÓXIMO bit
        bit.append(1 if (Vxp - Vxn) > 0 else 0)

    # 3. Bit de Transição C6
    v_step_c6 = (C6 / Ctot) * (Vref / 2)
    if bit[-1] == 1:
        Vxp -= v_step_c6
        Vxn += v_step_c6
    else:
        Vxp += v_step_c6
        Vxn -= v_step_c6
    bit.append(1 if (Vxp - Vxn) > 0 else 0)

    # 4. Array LSB (com atenuacao)
    # Nota: CtotalLSB e Cb devem estar acessíveis ou passados como argumento
    # Aqui vou usar a fórmula direta para facilitar
    atenuacao = Cb / (CtotalLSB + Cb)
    
    for j in range(Nlsb):
        v_step_lsb = (valores_CiLSB[j] / Ctot) * (Vref / 2) * atenuacao
        if bit[-1] == 1:
            Vxp -= v_step_lsb
            Vxn += v_step_lsb
        else:
            Vxp += v_step_lsb
            Vxn -= v_step_lsb
        bit.append(1 if (Vxp - Vxn) > 0 else 0)

    # Retornamos apenas os primeiros 12 bits para garantir o range 0-4095
    return bit[:12]

def converter_para_decimal(lista_bits):
    codigo_decimal = 0
    for b in lista_bits:
        codigo_decimal = codigo_decimal << 1
        codigo_decimal = codigo_decimal | b
    return codigo_decimal

for i in range(amostras):
    Vin = Vin_inicial[i]
    Vip = (Vref/2) + (Vin/2)
    Vin = (Vref/2) - (Vin/2)
    Vxp = Vip
    Vxn = Vin
    res_sar = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref)
    codigos.append(converter_para_decimal(res_sar))

plt.figure(figsize=(10,6))
plt.step(Vin_inicial, codigos)
plt.title("Característica de Transferência do ADC")
plt.xlabel("Vin (V)")
plt.ylabel("Código Decimal")
plt.grid(True)
plt.show()