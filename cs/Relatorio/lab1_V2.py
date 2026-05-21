import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Parametros do ADC
Nbits = 12
Vref = 1
Vlsb = (2 * Vref) / (2**Nbits)
codigos = []
codigos_real = []


multiplicadores_split_MSB = [32, 16, 8, 4, 2]
multiplicadores_split_C6 = [1]
multiplicadores_split_LSB = [16, 8, 4, 2, 1]

Nmsb = len(multiplicadores_split_MSB)
LC6 = len(multiplicadores_split_C6)
Nlsb = len(multiplicadores_split_LSB)


# --------------------------------------------------------------------------------------------------------------------------
# Funcao para gerar numeros aleatorios
def randn():
    return random.gauss(0, 1)


# Funcao para calcular os valores dos condensadores considerando os erros
def calcular_condensadores(SC, Cu, Cpb, Cpl, Cpm, Nlsb):
    valores_CiMSB = []
    for m in multiplicadores_split_MSB:
        unidades_a = (Cu/2) + (SC * np.random.randn(m))
        unidades_b = (Cu/2) + (SC * np.random.randn(m))
        valores_CiMSB.append(np.sum(unidades_a) + np.sum(unidades_b))

    valores_CiLSB = []
    for m in multiplicadores_split_LSB:
        unidades_lsb = Cu + (SC * np.random.randn(m))
        valores_CiLSB.append(np.sum(unidades_lsb))

    C6 = 1.0 * Cu

    # Para valores de Cb e Cdl:
    Cb = (2**Nlsb / ((2**Nlsb) - 1)) * Cu + Cpb
    #Cb = 2*Cu
    Cdl = (2**Nlsb - 1)*(Cb-Cu) - Cpl

    CtotalMSB = sum(valores_CiMSB) + Cpm
    CtotalLSB = sum(valores_CiLSB) + Cdl + Cpl

    C_ramo_LSB = (Cb * CtotalLSB) / (Cb + CtotalLSB)
    Ctot = CtotalMSB + C6 + C_ramo_LSB

    return valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB


# Funcao para simular o processo de conversao SAR
def sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl, v_off):
    bit = []
    vxp_hist = [Vxp]
    vxn_hist = [Vxn]

    # Conversao inicial
    if (Vxp-Vxn+v_off) < 0:
        bit.append(0)
    else:
        bit.append(1)

    # Conversao dos bits MSB
    for j in range(1, Nmsb+1):
        if bit[j-1] == 1:
            Vxp -= (valores_CiMSB[j-1] / Ctot) * Vref/2
            Vxn += (valores_CiMSB[j-1] / Ctot) * Vref/2
            if Vxp-Vxn+v_off > 0:
                bit.append(1)
            else:
                bit.append(0)
            
        elif bit[j-1] == 0:
            Vxp += (valores_CiMSB[j-1] / Ctot) * Vref/2
            Vxn -= (valores_CiMSB[j-1] / Ctot) * Vref/2

            if Vxp-Vxn+v_off > 0:
                bit.append(1)
            else:
                bit.append(0)

        vxp_hist.append(Vxp)
        vxn_hist.append(Vxn)

    # Conversao do bit C6
    if bit[5] == 1:
        Vxp -= (C6 / Ctot) * Vref/2
        Vxn += (C6 / Ctot) * Vref/2
        
        if Vxp-Vxn+v_off > 0:
            bit.append(1)
        else:
            bit.append(0)

    elif bit[5] == 0:
        Vxp += (C6 / Ctot) * Vref/2
        Vxn -= (C6 / Ctot) * Vref/2

        if Vxp-Vxn+v_off > 0:
            bit.append(1)
        else:
            bit.append(0)
    
    vxp_hist.append(Vxp)
    vxn_hist.append(Vxn)

    atenuacao = Cb / (CtotalLSB + Cpl + Cb)
    
    # Conversao dos bits LSB
    for j in range(7,Nlsb+7):
        if bit[j-1] == 1:
            Vxp -= (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao
            Vxn += (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao

            if Vxp-Vxn+v_off > 0:
                bit.append(1)
            else:
                bit.append(0)

        elif bit[j-1] == 0:
            Vxp += (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao
            Vxn -= (valores_CiLSB[j-7] / Ctot) * Vref/2 * atenuacao

            if Vxp-Vxn+v_off > 0:
                bit.append(1)
            else:
                bit.append(0)

        vxp_hist.append(Vxp)
        vxn_hist.append(Vxn)

    return bit[:12], vxp_hist, vxn_hist


# Funcao para converter a lista de bits em codigo decimal
def converter_para_decimal(lista_bits):
    codigo_decimal = 0
    for b in lista_bits:
        codigo_decimal = codigo_decimal << 1
        codigo_decimal = codigo_decimal | b
    return codigo_decimal


# Funcao calcular DNL e INL
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

    max_dnl_chip = np.max(np.abs(dnl))
    max_inl_chip = np.max(np.abs(inl))

    # Comentar para ver o Monte Carlo, senao aparecem muitos graficos

    # # Graficos
    # plt.figure(figsize=(10, 8))
    # plt.subplot(2,1,1)
    # plt.plot(codigos_apos[:-1], dnl)
    # plt.title("DNL por Codigo")
    # #plt.ylim(-0.5, 0.5)
    # plt.grid(True)
    
    # plt.subplot(2,1,2)
    # plt.plot(codigos_apos, inl)
    # plt.title("INL por Codigo")
    # #plt.ylim(-0.5, 0.5)
    # plt.grid(True)
    # plt.show()

    return  max_dnl_chip, max_inl_chip


# Funcao para calcular FFT e metricas de desempenho
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
    enob = (sndr - 1.76) / 6.02

    mag_db = 20 * np.log10(mag / np.max(mag))
    freq_norm = np.linspace(0, 0.5, len(mag_db))

    # Comentar para ver o Monte Carlo, senao aparecem muitos gráficos

    # # Grafico
    # if plot:
    #     # Tudo o que seja plt. deve estar indentado aqui dentro
    #     plt.figure(figsize=(10, 6))
    #     plt.plot(freq_norm, mag_db, color='blue', linewidth=0.7)
    #     plt.title(f"Dynamic Performance Analysis ({Nbits} bits)")
    #     plt.xlabel("Frequencia (MHz)")
    #     plt.ylabel("Magnitude (dB)")
    #     plt.ylim([-140, 5])
    #     plt.grid(True, which='both', linestyle='--', alpha=0.5)

    #     info_text = (f"Fin: {fin/1e6:.3f} MHz\n"
    #                 f"SNDR: {sndr:.2f} dB\n"
    #                 f"ENOB: {enob:.2f} bits")
    #     plt.text(0.65, 0.95, info_text, transform=plt.gca().transAxes,
    #             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7),
    #             fontsize=10, family='monospace')
    #     plt.show()
    # else:
    #     plt.close('all')

    return {"SNDR": sndr, "ENOB": enob}


# --------------------------------------------------------------------------------------------------------------------------
# Simulacao do processo de conversao para cada valor de Vin
num_codigos = 2**Nbits
pontos_por_lsb = 20
amostras = num_codigos * pontos_por_lsb

Vin_inicial = np.linspace(-Vref, Vref, amostras)

Cu = 1
SC = 0
Soffset = 0
Cpb = 0
Cpl = 0
Cpm = 0

valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(SC, Cu, Cpb, Cpl, Cpm, Nlsb)

for i in range(amostras):
    Vin = Vin_inicial[i]
    Vip = (Vref/2) + (Vin/2)
    Vin = (Vref/2) - (Vin/2)
    Vxp = Vip
    Vxn = Vin
    res_sar, h_vxp, h_vxn = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl, Soffset)
    codigos.append(converter_para_decimal(res_sar))

# Plot da caracteristica de transferencia do ADC
plt.figure(figsize=(10,6))
plt.step(Vin_inicial, codigos)
plt.title("Caracteristica de Transferencia do ADC")
plt.xlabel("Vin (V)")
plt.ylabel("Codigo Decimal")
plt.grid(True)
plt.show()

# --------------------------------------------------------------------------------------------------------------------------
# Simulacao do processo de conversao para cada valor de Vin mas com erro nos condensadores
SC_r = 0.03
Soffset_r = 0.01
Cpb_r = 1
Cpl_r = 1
Cpm_r = 1

valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(SC_r, Cu, Cpb_r, Cpl_r, Cpm_r, Nlsb)

for i in range(amostras):
    Vin = Vin_inicial[i]
    Vip = (Vref/2) + (Vin/2)
    Vin = (Vref/2) - (Vin/2)
    Vxp = Vip
    Vxn = Vin
    res_sar_real, h_vxp, h_vxn = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl, Soffset_r)
    codigos_real.append(converter_para_decimal(res_sar_real))

# Plot da caracteristica de transferencia do ADC
plt.figure(figsize=(10,6))
plt.step(Vin_inicial, codigos_real, color='red', label='ADC Real (com erro)', where='post')
plt.step(Vin_inicial, codigos, color='blue', label='ADC Ideal', where='post')
plt.title("Caracteristica de Transferencia do ADC")
plt.xlabel("Vin (V)")
plt.ylabel("Codigo Decimal")
plt.legend()
plt.grid(True)
plt.show()


# --------------------------------------------------------------------------------------------------------------------------
# Vxp e Vxn
Vin_exemplo = 0.1
Vip = (Vref/2) + (Vin_exemplo/2)
Vin = (Vref/2) - (Vin_exemplo/2)
Vxp = Vip
Vxn = Vin

res, h_vxp, h_vxn = sar_converter(Vip, Vin, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl, Soffset)

plt.figure(figsize=(12, 6))
plt.plot(h_vxp, 'o-', label='Vxp')
plt.plot(h_vxn, 'x-', label='Vxn')

posicoes = range(len(h_vxp))
labels_bits = [str(i+1) for i in posicoes]

plt.title(f"Evolucao de Vxp e Vxn para Vin = {Vin_exemplo}V")
plt.xlabel("Passo de Aproximacao (Bit)")
plt.ylabel("Tensao (V)")
plt.xticks(posicoes, labels=labels_bits)
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()


# --------------------------------------------------------------------------------------------------------------------------
# DNL e INL
calcular_dnl_inl(Vin_inicial, codigos, Nbits)

calcular_dnl_inl(Vin_inicial, codigos_real, Nbits)


# --------------------------------------------------------------------------------------------------------------------------
# Calculo de FFT e metricas de desempenho
N_fft = 2**(Nbits + 2)
fs = 100e6
M = 101
fin = (fs / N_fft) * M
t = np.arange(N_fft) / fs
Vin_seno = (Vref - Vlsb) * np.sin(2 * np.pi * fin * t)
codigos_fft = []

valores_CiMSB, valores_CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(SC, Cu, Cpb, Cpl, Cpm, Nlsb)

for v in Vin_seno:
    Vip = (Vref/2) + (v/2)
    Vin = (Vref/2) - (v/2)
    Vxp = Vip
    Vxn = Vin
    res_sar, h_vxp, h_vxn = sar_converter(Vxp, Vxn, Nmsb, Nlsb, valores_CiMSB, valores_CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl, Soffset)
    codigos_fft.append(converter_para_decimal(res_sar))

calcular_fft_e_metricas(codigos_fft, Nbits)


# --------------------------------------------------------------------------------------------------------------------------
# Monte Carlo
SC_mc = np.linspace(0, 0.0012, 10)
Numero_chips = 1000

sndr_results = np.zeros((len(SC_mc), Numero_chips))


for s_idx, sigma in enumerate(SC_mc):
    for i in range(Numero_chips):
        CiMSB, CiLSB, C6, Ctot, Cb, CtotalLSB = calcular_condensadores(sigma, Cu, Cpb, Cpl, Cpm, Nlsb)
        
        codigos_seno = []
        for v in Vin_seno:
            Vp = (Vref/2 + v/2)
            Vn = (Vref/2 - v/2)
            res, h_vxp, h_vxn = sar_converter(Vp, Vn, Nmsb, Nlsb, CiMSB, CiLSB, C6, Ctot, Vref, Cb, CtotalLSB, Cpl, Soffset)
            codigos_seno.append(converter_para_decimal(res))
        
        metrics = calcular_fft_e_metricas(codigos_seno, Nbits)
        sndr_results[s_idx, i] = metrics["SNDR"]
        
sigma_pct = SC_mc * 100

# Calculo das estatisticas finais
sndr_max = np.max(sndr_results, axis=1)
sndr_min = np.min(sndr_results, axis=1)
sndr_mean = np.mean(sndr_results, axis=1)

plt.figure(figsize=(12, 7))

# 1. Desenha a nuvem de pontos azuis
for s_idx, s_val in enumerate(sigma_pct):
    plt.scatter(np.ones(Numero_chips) * s_val, sndr_results[s_idx, :], 
                facecolors='none', edgecolors='tab:blue', alpha=0.3, s=10)

# 2. Desenha as linhas de tendencia
plt.plot(sigma_pct, sndr_max, color='blue', label='Max (Best Case)', linewidth=1.5)
plt.plot(sigma_pct, sndr_min, color='red', label='Min (Worst Case)', linewidth=1.5)
plt.plot(sigma_pct, sndr_mean, color='green', label='Mean', linewidth=2)

plt.title(r"SNDR vs. Capacitor Mismatch ($\sigma_{\Delta C/C}$)")
plt.xlabel(r"$\sigma$ of $\Delta C/C$ (%)")
plt.ylabel("SNDR (dB)")
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend()
plt.show()
 
# ── Parametros-alvo ──────────────────────────────────────────────────────────
ENOB_alvo  = 10.5      # bits
yield_alvo = 0.99      # 99 % dos chips devem cumprir o requisito
                       # → usamos o percentil 1 (pior 1 %) como criterio
 
# ── 1. ENOB a partir dos resultados Monte Carlo existentes ───────────────────
enob_results = (sndr_results - 1.76) / 6.02
 
enob_p1   = np.percentile(enob_results,  1, axis=1)   # critério yield 99 %
enob_mean = np.mean(enob_results,           axis=1)
enob_p99  = np.percentile(enob_results, 99, axis=1)
 
# ── 2. σ_max que satisfaz ENOB_alvo @ yield_alvo ────────────────────────────
try:
    f_interp  = interp1d(enob_p1, SC_mc, kind='linear', fill_value='extrapolate')
    sigma_max = float(f_interp(ENOB_alvo))
    sigma_max = np.clip(sigma_max, SC_mc[0], SC_mc[-1])
except Exception:
    sigma_max = SC_mc[np.argmin(np.abs(enob_p1 - ENOB_alvo))]
 
sigma_max_pct = sigma_max * 100      # em %
 
print(f"\n{'='*65}")
print(f"  Estimacao para ENOB ≥ {ENOB_alvo} bits  (yield ≥ {yield_alvo*100:.0f} %)")
print(f"{'='*65}")
print(f"  σ_max (ΔC/C)  = {sigma_max:.5f}  ({sigma_max_pct:.4f} %)")
 
# ── 3. area minima de Cu — modelo do paper (Fig. 4) ─────────────────────────
#
#   Recta da Fig. 4:   σ(%) = 0.408 / sqrt(Area_um2) + 0.01
#
#   Inverter para obter Area:
#       sigma_max_pct - 0.01 = 0.408 / sqrt(Area)
#       sqrt(Area) = 0.408 / (sigma_max_pct - 0.01)
#       Area = (0.408 / (sigma_max_pct - 0.01))^2
#
#   Se sigma_max_pct <= 0.01 (abaixo do patamar do modelo), area seria infinita
#   — nesse caso usa-se o maior valor do eixo x da Fig. 4 (area ≈ 16 um^2).
 
a_fit   = 0.408   # coeficiente da recta y = a/sqrt(Area) + b  (Fig. 4)
b_fit   = 0.01    # offset
 
delta = sigma_max_pct - b_fit
if delta > 1e-6:
    Area_Cu_um2 = (a_fit / delta) ** 2
else:
    Area_Cu_um2 = (a_fit / (0.112 - b_fit)) ** 2   # fallback: ponto (16, 0.112)
    print(f"  AVISO: σ_max muito baixo; usando Area_Cu = 16 um^2 como limite superior.")
 
# Lado do condensador quadrado equivalente
lado_Cu = np.sqrt(Area_Cu_um2)
 
print(f"\n  Modelo mismatch (Fig. 4):  σ(%) = {a_fit}/√Area + {b_fit}")
print(f"  Area minima de Cu         = {Area_Cu_um2:.2f} um^2  ({lado_Cu:.2f}×{lado_Cu:.2f} um^2)")
 
# Cu real do paper para comparacao
Area_Cu_paper = 4 * 4          # um^2
sigma_paper   = a_fit / np.sqrt(Area_Cu_paper) + b_fit
print(f"\n  [Referencia paper]  Cu = 4×4 um^2  →  σ = {sigma_paper:.3f} %")
 
# ── 4. Capacidade de Cu e total do DAC ──────────────────────────────────────
#
#   O paper indica:
#     - Cu_paper = 6.5 fF  para  Area = 16 um^2
#     - Densidade: C_ox = 6.5 / 16 = 0.406 fF/um^2
#
#   Escalar linearmente para a area estimada.
 
C_ox_fF_um2 = 6.5 / 16.0          # fF/um^2  (do paper, texto Sec. III-B)
Cu_fF        = Area_Cu_um2 * C_ox_fF_um2
 
# Número total de unidades Cu numa array single-ended (da Tabela 1 do paper):
#   MSB (bits 1-7, com redundancia):  56+32+16+10+6+4+2 = 126 Cu
#   C6 (bridge):  ≈ 2 Cu  (CB = 2Cu conforme paper)
#   LSB (bits 8-14, com redundancia): 1+7+4+2+1+1+15 = 31 Cu  (aprox., ver Tabela 1)
#   C_total single-ended paper ≈ 840 fF / 6.5 fF ≈ 129 Cu  → usa valor do paper
 
N_Cu_single   = round(840 / 6.5)   # ≈ 129 unidades Cu por ramo single-ended
N_Cu_total    = 2 * N_Cu_single     # diferencial (p + n)
 
C_total_fF    = N_Cu_total * Cu_fF
 
# Escalar a área total do DAC
Area_DAC_um2  = N_Cu_total * Area_Cu_um2
 
print(f"\n  Densidade C_ox            = {C_ox_fF_um2:.3f} fF/µm²")
print(f"  Cu estimado               = {Cu_fF:.2f} fF")
print(f"  N_Cu (single-ended)       ≈ {N_Cu_single}  unidades Cu")
print(f"  C_total DAC (diferencial) ≈ {C_total_fF:.0f} fF")
print(f"  Area total DAC            ≈ {Area_DAC_um2:.0f} µm²  ({Area_DAC_um2*1e-6:.4f} mm²)")
 
# ── 5. Energia por conversao (DAC) ──────────────────────────────────────────
#
#   Para um SAR com switching monotonico/hibrido, a energia media de
#   chaveamento do DAC e aproximadamente:
#       E_DAC ≈ (1/2) · C_total · Vref^2    (limite superior conservador)
#
#   Vref = 1.2 V  (paper, Sec. II)
#
#   FoM_DAC = E_DAC / 2^ENOB
 
Vref_paper   = 1.2       # V
alpha        = 0.5       # factor de chaveamento medio SAR
E_DAC_fJ     = alpha * C_total_fF * (Vref_paper ** 2)
FoM_DAC      = E_DAC_fJ / (2 ** ENOB_alvo)
 
# Comparacao com o paper (ENOB=10.17, P=0.8mW, fs=100MHz)
E_total_paper_fJ = (0.8e-3 / 100e6) * 1e15     # fJ por conversao (total)
FoM_paper        = E_total_paper_fJ / (2**10.17)
 
print(f"\n  Vref                      = {Vref_paper} V")
print(f"  Energia DAC estimada      ≈ {E_DAC_fJ:.1f} fJ/conversao")
print(f"  FoM DAC estimado          ≈ {FoM_DAC:.2f} fJ/conversion-step")
print(f"\n  [Referencia paper]  E_total = {E_total_paper_fJ:.1f} fJ  →  FoM = {FoM_paper:.2f} fJ/conv-step")
print(f"{'='*65}\n")
 
# ── 6. Gráfico ENOB vs σ com marcação do ponto-alvo ─────────────────────────
sigma_pct_plot = SC_mc * 100
 
fig, ax = plt.subplots(figsize=(10, 6))
 
for s_idx, s_val in enumerate(sigma_pct_plot):
    ax.scatter(np.ones(Numero_chips) * s_val, enob_results[s_idx, :],
               facecolors='none', edgecolors='tab:blue', alpha=0.2, s=8)
 
ax.plot(sigma_pct_plot, enob_p99,  color='blue',  label='Percentil 99 (melhor)',          linewidth=1.5)
ax.plot(sigma_pct_plot, enob_mean, color='green', label='Média',                           linewidth=2)
ax.plot(sigma_pct_plot, enob_p1,   color='red',   label='Percentil 1  (critério yield 99%)', linewidth=1.5)
 
ax.axhline(ENOB_alvo,       color='black',      linestyle='--', linewidth=1.2,
           label=f'Alvo: ENOB = {ENOB_alvo} bits')
ax.axvline(sigma_max_pct,   color='darkorange', linestyle=':',  linewidth=1.5,
           label=f'σ_max = {sigma_max_pct:.3f} %  (yield {yield_alvo*100:.0f} %)')
 
# Ponto do paper
ax.axvline(sigma_paper, color='purple', linestyle='-.', linewidth=1.2,
           label=f'Paper: Cu 4×4 µm²  →  σ = {sigma_paper:.3f} %')
 
ax.set_title(r"ENOB vs. Capacitor Mismatch ($\sigma_{\Delta C/C}$) — Alínea h)")
ax.set_xlabel(r"$\sigma$ de $\Delta C/C$ (%)")
ax.set_ylabel("ENOB (bits)")
ax.legend(fontsize=9)
ax.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()