# Aba 5 — Engenharia, Formulação Matemática e Rastreabilidade de Parâmetros
# Documentação técnica sincronizada com os parâmetros e calibrações reais de modelo.jl.

function cartao_engenharia!(gl, titulo, itens, cor_titulo, largura_texto::Int)
    NL = string('\n')
    try
        Box(gl[1:2, 1]; color = CARTAO, strokecolor = BORDA, strokewidth = 1.4, cornerradius = 10)
    catch
        Box(gl[1:2, 1]; color = CARTAO, strokecolor = BORDA, strokewidth = 1.4)
    end
    Label(gl[1, 1], titulo;
          fontsize = 14, font = :bold, color = cor_titulo, halign = :left,
          padding = (16, 16, 10, 4), word_wrap = true, width = largura_texto)
    Label(gl[2, 1], join(itens, NL);
          fontsize = 12, color = TEXTO, halign = :left, valign = :top, justification = :left,
          padding = (16, 16, 2, 10), word_wrap = true, width = largura_texto)
    rowsize!(gl, 1, Auto())
    rowsize!(gl, 2, Auto())
    rowgap!(gl, 4)
    return
end

function linhas_rastreabilidade(d, filtro::Int)
    rows = [
        (1, "PAM alvo (ATHOS-3)", f0(PAM_ALVO_ATHOS), "mmHg", "Critério clínico · Khanna et al. (2017)"),
        (1, "PAM referência (SSC)", f0(PAM_ALVO_SSC), "mmHg", "Diretriz sepse · Evans et al. (2021)"),
        (1, "Pressão venosa (Pᵥₑₙ)", f0(P_VEN), "mmHg", "Referência basal · Melo & Peters (1999)"),
        (1, "Impedância aórtica (Z꜀)", "0,06 × Rₚ", "mmHg·s/ml", "Relação proximal · Melo & Peters (1999)"),
        (1, "Complacência adulto (C₀)", f2(C_ADULTO), "ml/mmHg", "Planta basal · Melo & Peters (1999)"),
        (1, "Fator Rₚ idade (70 a)", f2(PARAM.frp[]), "×", "Enrijecimento · Ferrari et al. (2003)"),
        (1, "Fator C idade (70 a)", f2(PARAM.fcx[]), "×", "Perda elástica · Ferrari et al. (2003)"),
        (1, "FC no choque séptico", f0(PARAM.fc[]), "bpm", "Taquicardia · Khanna et al. (2017)"),
        (1, "DC choque idoso", f1(PARAM.dci[]), "L/min", "Calibração · ATHOS-3 / RVS 455"),
        (2, "Eₘₐₓ Angiotensina II", f2(PARAM.Emax[]), "adim.", "Sensibilidade · Leisman et al. (2024)"),
        (2, "EC₅₀ Angiotensina II", f1(d.EC50), "ng/kg/min", "Calibrado p/ Eₘₐₓ · Leisman et al. (2024)"),
        (2, "Fator r do BRA (AT₁)", f1(d.rb), "× EC₅₀", "Subgrupo n = 22 · Leisman et al. (2024)"),
        (2, "Fator r do IECA", "1,0", "× EC₅₀", "Subgrupo n = 29 · Leisman et al. (2024)"),
        (2, "Meia-vida Ang II (t₁/₂)", f1(T_MEIA_AII), "min", "Farmacocinética · Basso & Terragno"),
        (2, "Dose inicial / máx Ang II", f0(DOSE_INICIAL) * " / " * f0(DOSE_MAX), "ng/kg/min", "Protocolo ATHOS-3 · Khanna et al. (2017)"),
        (2, "Cortes prognósticos F₁₅", f2(F15_CORTES[1]) * " – " * f2(F15_CORTES[2]), "adim.", "Déficit em 15 min · Leisman et al. (2024)"),
        (3, "Ganho de titulação", f1(PARAM.ganho[]), "ng/kg/min/mmHg", "Passo proporcional · Protocolo grupo"),
        (3, "Horizonte de simulação", f0(HORIZONTE), "min", "Janela clínica · 3 h (ATHOS-3)"),
        (3, "Passo Euler batimento", replace(@sprintf("%.1f", DT_EULER * 1000), "." => ","), "ms", "Integração EDO · Estabilidade numérica"),
        (3, "Eₘₐₓ relativo NE", f2(PARAM.ne_frac[]), "× Eₘₐₓ,AII", "Premissa do grupo (não calibrada)"),
        (3, "EC₅₀ Noradrenalina", f2(PARAM.ne_ec50[]), "µg/kg/min", "Premissa do grupo (dose alta NE)"),
        (3, "Meia-vida NE (t₁/₂)", f1(PARAM.ne_t12[]), "min", "Premissa do grupo (desmame :troca)")
    ]

    filtradas = filtro == 0 ? rows : [r for r in rows if r[1] == filtro]
    cab = rpad("Parâmetro", 28) * lpad("Valor", 12) * "   " * rpad("Unidade", 16) * "Origem / Referência"
    sep = repeat("─", 94)
    corpo_linhas = [rpad(r[2], 28) * lpad(r[3], 12) * "   " * rpad(r[4], 16) * r[5] for r in filtradas]
    return join(vcat([cab, sep], corpo_linhas), string('\n'))
end

function aba_engenharia(corpo)
    corpo.valign = :top
    FILTRO_TAB = Observable(0)
    fonte_mono = fonte_ou("Consolas", :regular)
    larg_esq = 730
    larg_txt = larg_esq - 32

    # ================================================================
    # COLUNA ESQUERDA: EQUAÇÕES E MÉTODOS NUMÉRICOS
    # ================================================================
    esq = corpo[1, 1] = GridLayout(valign = :top)

    c1 = esq[1, 1] = GridLayout()
    cartao_engenharia!(c1,
        "1 · Planta Hemodinâmica — Windkessel de 3 Elementos (WK3)",
        [
            "• Conservação de massa no reservatório arterial:   C · (dP꜀/dt) = Qᵢₙ(t) − (P꜀ − Pᵥₑₙ) / Rₚ",
            "• Pressão aórtica instantânea:   Pₐₒ(t) = P꜀(t) + Z꜀ · Qᵢₙ(t),   com Z꜀ = 0,06 · Rₚ",
            "• Ejeção ventricular sistólica:   Qᵢₙ(t) = (π·VS / 2·Tₑⱼ) · sin(π·t / Tₑⱼ)   para 0 ≤ t < Tₑⱼ",
            "• Envelhecimento vascular (70 anos):   Rₚ,ᵢ = f_Rp · Rₚ,₀   e   Cᵢ = f_C · C₀ (↑ pressão de pulso)"
        ],
        AZUL, larg_txt)

    c2 = esq[2, 1] = GridLayout()
    cartao_engenharia!(c2,
        "2 · Farmacocinética (PK), Farmacodinâmica (Hill) & Controle",
        [
            "• Compartimento PK (1ª ordem):   dx/dt = kₑ · (u(t) − x),   com kₑ = ln(2) / t₁/₂",
            "• Ação vasoconstritora (Hill):   E_AII(x) = Eₘₐₓ · x / (r · EC₅₀ + x)   (r = 1 sem BRA; r > 1 com BRA)",
            "• Modulação da resistência periférica:   Rₚ(t) = Rₚ,₀ · (1 + E_NE(t) + E_AII(t))",
            "• Titulação ATHOS-3: ajuste discreto a cada 5 min (0–15 min) e 15 min (15–180 min) até PAM ≥ 75 mmHg"
        ],
        LARANJA, larg_txt)

    c3 = esq[3, 1] = GridLayout()
    cartao_engenharia!(c3,
        "3 · Métodos Numéricos, Aproximações e Escopo de Engenharia",
        [
            "• Integração intra-ciclo: método de Euler explícito com passo fixo dt = 0,1 ms (10.000 Hz)",
            "• Horizonte de 3 horas: aproximação quase-estática PAM(t) = (Rₚ(t) + Z꜀(t))·DC + Pᵥₑₙ (passo 1 s)",
            "• Desmame da noradrenalina (:troca): premissa exploratória do grupo (redução gradual após PAM ≥ 75)",
            "• Pendência futura: controlador PI contínuo com split-range coordenado entre Ang II e Noradrenalina"
        ],
        VERDE, larg_txt)

    rowgap!(esq, 12)

    # ================================================================
    # COLUNA DIREITA: TABELA DE RASTREABILIDADE (MOLDURA ABRAÇA O TEXTO)
    # ================================================================
    dir = corpo[1, 2] = GridLayout(valign = :top)
    try
        Box(dir[1:4, 1]; color = PAINEL, strokecolor = BORDA, strokewidth = 1.4, cornerradius = 10)
    catch
        Box(dir[1:4, 1]; color = PAINEL, strokecolor = BORDA, strokewidth = 1.4)
    end

    Label(dir[1, 1], "4 · Tabela de Rastreabilidade dos Parâmetros (tempo real)";
          fontsize = 14, font = :bold, color = AMARELO, halign = :left,
          padding = (18, 18, 10, 2), tellwidth = false)

    filtros_gl = dir[2, 1] = GridLayout(halign = :left, padding = (18, 18, 2, 4))
    nomes_f = ["Todos (22)", "1 · Windkessel (WK3)", "2 · Ang II & BRA", "3 · Controle & NE"]
    bts_f = Any[]
    for (idx, rot) in enumerate(nomes_f)
        val = idx - 1
        bt = Button(filtros_gl[1, idx];
                    label = rot,
                    width = 175,
                    fontsize = 11,
                    buttoncolor = FILTRO_TAB[] == val ? ATIVA : INATIVA,
                    labelcolor = FILTRO_TAB[] == val ? :black : :white)
        on(_ -> (FILTRO_TAB[] = val), bt.clicks)
        push!(bts_f, bt)
    end
    on(FILTRO_TAB) do f
        for (idx, bt) in enumerate(bts_f)
            ativo = (idx - 1) == f
            bt.buttoncolor[] = ativo ? ATIVA : INATIVA
            bt.labelcolor[]  = ativo ? :black : :white
        end
    end
    colgap!(filtros_gl, 8)

    txt_tabela = lift((d, f) -> linhas_rastreabilidade(d, f), DADOS, FILTRO_TAB)
    Label(dir[3, 1], txt_tabela;
          font = fonte_mono, fontsize = 11, color = TEXTO, halign = :left, valign = :top,
          justification = :left, padding = (18, 18, 4, 4), tellwidth = false)

    Label(dir[4, 1], "Nota: os valores numéricos acima refletem automaticamente qualquer ajuste feito nos sliders das abas 1, 2, 3 e 4.";
          fontsize = 11, color = GRIS, halign = :left, padding = (18, 18, 2, 10),
          tellwidth = false, word_wrap = true)

    rowsize!(dir, 1, Auto())
    rowsize!(dir, 2, Auto())
    rowsize!(dir, 3, Auto())
    rowsize!(dir, 4, Auto())
    rowgap!(dir, 4)

    colsize!(corpo, 1, Fixed(larg_esq))
    colgap!(corpo, 16)
    return
end