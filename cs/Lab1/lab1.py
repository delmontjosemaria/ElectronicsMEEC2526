import math
import random
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft


# Parâmetros do ADC
Nbits = 12
Vref = 1.0
Vlsb = (2 * Vref) / (2**Nbits)
codigos = []


# Parâmetros dos condensadores
SC = 0e-2
Soffset = 0e-3
Cu = 1.0
Cpb = 0
Cpl = 0
Cpm = 0

multiplicadores_split_MSB = [32, 16, 8, 4, 2]
multiplicadores_split_C6 = [1]
multiplicadores_split_LSB = [16, 8, 4, 2, 1]

Nmsb = len(multiplicadores_split_MSB)
LC6 = len(multiplicadores_split_C6)
Nlsb = len(multiplicadores_split_LSB)

def randn():
    return random.gauss(0, 1)

def calcular_condensadores(SC, Cu, Cpb, Cpl, Cpm, Nlsb):
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

    return valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB


# Simulação do processo de conversão
amostras = 10000
Vin_inicial = np.linspace(-Vref, Vref, amostras)

# Função para simular o processo de conversão SAR
def sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl):
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


# Simulação do processo de conversão para cada valor de Vin
valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(SC, Cu, Cpb, Cpl, Cpm, Nlsb)

for i in range(amostras):
    Vin = Vin_inicial[i]
    Vip = (Vref/2) + (Vin/2)
    Vin = (Vref/2) - (Vin/2)
    Vxp = Vip
    Vxn = Vin
    res_sar = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl)
    codigos.append(converter_para_decimal(res_sar))


# Plot da característica de transferência do ADC
plt.figure(figsize=(10,6))
plt.step(Vin_inicial, codigos)
plt.title("Característica de Transferência do ADC")
plt.xlabel("Vin (V)")
plt.ylabel("Código Decimal")
plt.grid(True)
plt.show()

# --------------------------------------------------------------------------------------------------------------------------
# DNL e INL
def calcular_dnl_inl(vin, codigos, nbits):
    codigo = np.array(codigos)
    vin = np.array(vin)

    indices_mudanca = np.where(np.diff(codigo) > 0)[0]
    v_transicoes = vin[indices_mudanca]
    codigos_apos = codigo[indices_mudanca + 1]

    v_zero = v_transicoes[0]
    v_max = v_transicoes[-1]
    num_degraus = codigos_apos[-1] - codigos_apos[0]

    vlsb_real = (v_max - v_zero) / num_degraus

    v_ideais = v_zero + (codigos_apos - codigos_apos[0]) * vlsb_real
    larguras_reais = np.diff(v_transicoes)

    inl = (v_transicoes - v_ideais) / vlsb_real
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

    return

calcular_dnl_inl(Vin_inicial, codigos, Nbits)

# --------------------------------------------------------------------------------------------------------------------------
# Função para calcular FFT e métricas de desempenho
def calcular_fft_e_metricas(codigos, Nbits, plot=True):
    N = len(codigos)
    
    janela = np.blackman(N)
    sinal_ac = (np.array(codigos) - np.mean(codigos)) * janela

    espetro = np.fft.fft(sinal_ac)
    mag = np.abs(espetro)[:N//2]
    psd = mag**2

    sinal = np.argmax(psd)
    p_s = np.sum(psd[max(0, sinal-2) : min(N//2, sinal+3)])
    p_d = 0

    h_max = int((N/2) / sinal) if sinal > 0 else 0
    
    for i in range(2, h_max + 1):
        h_idx_teorica = i * sinal
        fold = h_idx_teorica % N

        if fold >= N // 2:
            h_idx_real = N - fold
        else:
            h_idx_real = fold
            
        h_idx_final = int(h_idx_real)

        if h_idx_final > 2 and h_idx_final < (N // 2 - 3):
            p_d += np.sum(psd[h_idx_final-2 : h_idx_final+3])

    p_total = np.sum(psd)
    p_n = p_total - p_s - p_d
    p_n = max(p_n, 1e-20)

    sndr = 10 * np.log10(p_s / (p_n + p_d))
    thd = 10 * np.log10(p_d / p_s) if p_d > 0 else -140
    enob = (sndr - 1.76) / 6.02

    psd_sem_sinal = np.copy(psd)
    psd_sem_sinal[max(0, sinal-5) : min(N//2, sinal+6)] = 0
    p_spur_max = np.max(psd_sem_sinal)
    sfdr = 10 * np.log10(p_s / p_spur_max)

    mag_db = 20 * np.log10(mag / np.max(mag))
    freq_norm = np.linspace(0, 0.5, len(mag_db))

    # Gráfico
    if plot:
        # Tudo o que seja plt. deve estar indentado aqui dentro
        plt.figure(figsize=(10, 6))
        plt.plot(freq_norm, mag_db, color='blue', linewidth=0.7)
        plt.title(f"Dynamic Performance Analysis ({Nbits} bits)")
        plt.xlabel("Frequência Normalizada (f / fs)")
        plt.ylabel("Magnitude (dB)")
        plt.ylim([-140, 5])
        plt.grid(True, which='both', linestyle='--', alpha=0.5)

        info_text = (f"SNDR: {sndr:.2f} dB\n"
                    f"ENOB: {enob:.2f} bits\n"
                    f"THD:  {thd:.2f} dB\n"
                    f"SFDR: {sfdr:.2f} dB")
        plt.text(0.32, -30, info_text, bbox=dict(facecolor='white', alpha=0.8), fontsize=10, family='monospace')
        plt.show()
    else:
        plt.close('all')

    return {"SNDR": sndr, "ENOB": enob, "THD": thd, "SFDR": sfdr}

# --------------------------------------------------------------------------------------------------------------------------
# Cálculo de FFT e métricas de desempenho
N_fft = 4096
M = 13
fs = 1000
fin = (M / N_fft) * fs
t = np.arange(N_fft) / fs
Vin_seno = (Vref - Vlsb) * np.sin(2 * np.pi * fin * t)
codigos_fft = []

valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(SC, Cu, Cpb, Cpl, Cpm, Nlsb)

for v in Vin_seno:
    Vip = (Vref/2) + (v/2)
    Vin = (Vref/2) - (v/2)
    Vxp = Vip
    Vxn = Vin
    res_sar = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl)
    codigos_fft.append(converter_para_decimal(res_sar))

calcular_fft_e_metricas(codigos_fft, Nbits)

# --------------------------------------------------------------------------------------------------------------------------
# Monte Carlo
num_chips = 50
resultados_enob = []
SC_MC = 0.001

for i in range(num_chips):
    valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(SC_MC, Cu, Cpb, Cpl, Cpm, Nlsb)

    codigos_mc = []
    for v in Vin_seno:
        Vp = (Vref/2 + v/2)
        Vn = (Vref/2 - v/2)
        res = sar_converter(Vp, Vn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl)
        codigos_mc.append(converter_para_decimal(res))
    
    metrics = calcular_fft_e_metricas(codigos_mc, Nbits, plot=False)
    resultados_enob.append(metrics["ENOB"])

plt.figure(figsize=(9, 6))
plt.hist(resultados_enob, bins=15, color='skyblue', edgecolor='black')
plt.axvline(np.mean(resultados_enob), color='red', linestyle='dashed', label=f'Média: {np.mean(resultados_enob):.2f}')
plt.title(f"Análise de Monte Carlo")
plt.xlabel("ENOB (bits)")
plt.ylabel("Frequência")
plt.legend()
plt.show()
