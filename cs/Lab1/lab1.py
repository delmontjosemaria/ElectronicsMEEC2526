"""import math
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

weightsMSB = [32, 16, 8, 4, 2]
weightC6 = [1]
weightsLSB = [16, 8, 4, 2, 1]

nMSB = len(weightsMSB)
nLSB = len(weightsLSB)

def randn():
    return random.gauss(0, 1)

valueCia = [m * (Cu + SC * randn())/2 for m in weightsMSB]
valueCib = [m * (Cu + SC * randn())/2 for m in weightsMSB]
valueCiMSB = [valueCia[i] + valueCib[i] for i in range(len(valueCia))]
valueCiLSB = [m * (Cu + SC * randn()) for m in weightsLSB]

C6 = 1.0 * Cu
Cb = (2*nLSB / (2*nLSB - 1)) * Cu + Cpb
Cdl = (2*nLSB - 1)*(Cb-Cu) - Cpl

CtotalMSB = sum(valueCia) + sum(valueCib) + Cpm
CtotalLSB = sum(valueCiLSB) + Cdl + Cpl

C_ramo_LSB = (Cb * (CtotalLSB)) / (Cb + CtotalLSB)
Ctot = CtotalMSB + C6 + C_ramo_LSB


# Simulação do processo de conversão
amostras = 8000
Vin_inicial = np.linspace(-Vref, Vref, amostras)

def comparator(vxp:float, vxn:float):
    return 1 if vxp > vxn else 0

def sarConverter(Vxp, Vxn, nMSB, nLSB, valueCiMSB, valueCiLSB, C6, Ctot, Vref):
    bit = []
    
    bit.append(comparator(Vxp, Vxn))

    for j in range(nMSB):
        # Usamos o ÚLTIMO bit decidido para ajustar a tensão
        vStepMSB = (valueCiMSB[j] / Ctot) * (Vref / 2)
        if bit[-1] == 1:
            Vxp -= vStepMSB
            Vxn += vStepMSB
        else:
            Vxp += vStepMSB
            Vxn -= vStepMSB
        
        # Agora decidimos o PRÓXIMO bit
        bit.append(comparator(Vxp, Vxn))

    # 3. Bit de Transição C6
    vStepC6 = (C6 / Ctot) * (Vref / 2)
    if bit[-1] == 1:
        Vxp -= vStepC6
        Vxn += vStepC6
    else:
        Vxp += vStepC6
        Vxn -= vStepC6
    bit.append(comparator(Vxp, Vxn))

    # 4. Array LSB (com atenuacao)
    # Nota: CtotalLSB e Cb devem estar acessíveis ou passados como argumento
    # Aqui vou usar a fórmula direta para facilitar
    atenuacao = Cb / (CtotalLSB + Cb)
    
    for j in range(nLSB):
        vStepLSB = (valueCiLSB[j] / Ctot) * (Vref / 2) * atenuacao
        if bit[-1] == 1:
            Vxp -= vStepLSB
            Vxn += vStepLSB
        else:
            Vxp += vStepLSB
            Vxn -= vStepLSB
        bit.append(comparator(Vxp, Vxn))

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
    res_sar = sarConverter(Vxp, Vxn, nMSB, nLSB, valueCiMSB, valueCiLSB, C6, Ctot, Vref)
    codigos.append(converter_para_decimal(res_sar))

plt.figure(figsize=(10,6))
plt.step(Vin_inicial, codigos)
plt.title("Característica de Transferência do ADC")
plt.xlabel("Vin (V)")
plt.ylabel("Código Decimal")
plt.grid(True)
plt.show()"""

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

#multiplicadores_split_MSB = [16, 8, 4, 2, 1]
multiplicadores_split_MSB = [32, 16, 8, 4, 2]
multiplicadores_split_C6 = [1]
multiplicadores_split_LSB = [16, 8, 4, 2, 1]

Nmsb = len(multiplicadores_split_MSB)
LC6 = len(multiplicadores_split_C6)
Nlsb = len(multiplicadores_split_LSB)

def randn():
    return random.gauss(0, 1)

valores_Cia = [m * (Cu + SC * randn())/2 for m in multiplicadores_split_MSB]
valores_Cib = [m * (Cu + SC * randn())/2 for m in multiplicadores_split_MSB]
valores_CiMSB = [valores_Cia[i] + valores_Cib[i] for i in range(len(valores_Cia))]
valores_CiLSB = [m * (Cu + SC * randn()) for m in multiplicadores_split_LSB]

C6 = 1.0 * Cu
Cb = (2**Nlsb / ((2**Nlsb) - 1)) * Cu + Cpb
Cdl = (2**Nlsb - 1)*(Cb-Cu) - Cpl

CtotalMSB = sum(valores_Cia) + sum(valores_Cib) + Cpm
CtotalLSB = sum(valores_CiLSB) + Cdl + Cpl

C_ramo_LSB = (Cb * CtotalLSB) / (Cb + CtotalLSB)
Ctot = CtotalMSB + C6 + C_ramo_LSB


# Simulação do processo de conversão
amostras = 10000
Vin_inicial = np.linspace(-Vref, Vref, amostras)

# Função para simular o processo de conversão SAR
def sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref):
    bit = []

    # Conversão inicial
    if (Vxp - Vxn) < 0:
        bit.append(0)
    else:
        bit.append(1)

    # Conversão dos bits MSB
    for j in range(1, Nmsb+1):
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

    # Conversão do bit C6
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
    
    # Conversão dos bits LSB
    for j in range(7,Nlsb+7):
        if bit[j-1] == 1:
            Vxp -= (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao
            Vxn += (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao

            if Vxp-Vxn > 0:
                bit.append(1)
            else:
                bit.append(0)

        elif bit[j-1] == 0:
            Vxp += (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao
            Vxn -= (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao

            if Vxp-Vxn > 0:
                bit.append(1)
            else:
                bit.append(0)


    return bit[:12]

# Função para converter a lista de bits em código decimal
def converter_para_decimal(lista_bits):
    codigo_decimal = 0
    for b in lista_bits:
        codigo_decimal = codigo_decimal << 1
        codigo_decimal = codigo_decimal | b
    return codigo_decimal

# DNL e INL
def calcular_dnl_inl(vin, codigos, nbits):
    codigo = np.array(codigos)
    vin = np.array(vin)

    # Encontrar os índices onde o código muda (transições)
    indices_mudanca = np.where(np.diff(codigo) > 0)[0]

    # Tensões onde ocorreram as transições
    v_transicoes = vin[indices_mudanca]
    
    # Códigos correspondentes após a transição
    codigos_apos = codigo[indices_mudanca + 1]

    # --- Cálculo de Vzero e Vmax ---
    v_zero = v_transicoes[0]
    v_max = v_transicoes[-1]

    num_degraus = codigos_apos[-1] - codigos_apos[0]

    # --- VLSB REAL ---
    vlsb_real = (v_max - v_zero) / num_degraus

    # --- DNL e INL ---
    v_ideais = v_zero + (codigos_apos - codigos_apos[0]) * vlsb_real
    inl = (v_transicoes - v_ideais) / vlsb_real
    
    larguras_reais = np.diff(v_transicoes)
    dnl = (larguras_reais / vlsb_real) - 1
    
    # Gráficos
    plt.figure(figsize=(10, 8))
    plt.subplot(2,1,1)
    plt.plot(codigos_apos[:-1], dnl)
    plt.title("DNL por Código")
    plt.grid(True)
    
    plt.subplot(2,1,2)
    plt.plot(codigos_apos, inl)
    plt.title("INL por Código")
    plt.grid(True)
    plt.show()

    return dnl, inl


# Simulação do processo de conversão para cada valor de Vin
for i in range(amostras):
    Vin = Vin_inicial[i]
    Vip = (Vref/2) + (Vin/2)
    Vin = (Vref/2) - (Vin/2)
    Vxp = Vip
    Vxn = Vin
    res_sar = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref)
    codigos.append(converter_para_decimal(res_sar))


# Plot da característica de transferência do ADC
plt.figure(figsize=(10,6))
plt.step(Vin_inicial, codigos)
plt.title("Característica de Transferência do ADC")
plt.xlabel("Vin (V)")
plt.ylabel("Código Decimal")
plt.grid(True)
plt.show()

calcular_dnl_inl(Vin_inicial, codigos, Nbits)