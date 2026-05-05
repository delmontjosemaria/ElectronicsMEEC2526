def conversor_sar_monotonico(v_ip, v_in, v_ref_p, v_ref_n, resolucao=12):
    """
    Simula um ADC SAR com chaveamento monotônico.
    """
    # 1. Definir a tensão de referência total e o tamanho do primeiro passo (MSB)
    v_ref_total = v_ref_p - v_ref_n
    
    # No chaveamento monotônico, o primeiro degrau afeta a tensão em Vref/2
    passo_tensao = v_ref_total / 2.0 
    
    # 2. Amostragem (Fase 1 do circuito)
    v_xp = v_ip
    v_xn = v_in
    
    # Array para guardar os nossos bits
    codigo_digital = []
    
    # 3. O Loop Assíncrono / Motor de Busca
    for i in range(resolucao):
        
        # Comparador avalia quem é maior
        if v_xp > v_xn:
            # Decisão é 1
            codigo_digital.append(1)
            # Como a rede P era maior, a chave do DAC P muda, baixando V_xp
            v_xp = v_xp - passo_tensao
        else:
            # Decisão é 0
            codigo_digital.append(0)
            # Como a rede N era maior (ou igual), a chave do DAC N muda, baixando V_xn
            v_xn = v_xn - passo_tensao
            
        # Prepara o DAC para o próximo bit (próximo capacitor é metade do tamanho)
        passo_tensao = passo_tensao / 2.0
        
    return codigo_digital

# ==========================================
# Testando o nosso ADC Virtual
# ==========================================
# Tensões de Referência (ex: 3.3V e 0V)
V_REFP = 3.3
V_REFN = 0.0

# Tensões de Entrada (O nosso "x")
# Vamos testar uma entrada onde V_ip é ligeiramente maior que V_in
V_IP_teste = 2.1
V_IN_teste = 1.4 

# Executar a conversão para obter o "y"
bits_y = conversor_sar_monotonico(V_IP_teste, V_IN_teste, V_REFP, V_REFN)

print(f"Tensão Diferencial de Entrada: {V_IP_teste - V_IN_teste:.4f} V")
print(f"Código Digital Resultante (12-bits): {bits_y}")

# Opcional: Converter array binário para inteiro (valor de 0 a 4095)
valor_inteiro = sum(bit * (2 ** (11 - i)) for i, bit in enumerate(bits_y))
print(f"Valor Decimal: {valor_inteiro} de 4095")