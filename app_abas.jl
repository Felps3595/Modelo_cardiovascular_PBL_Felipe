# Plataforma em abas do choque vasodilatador refratário (GLMakie).
#   Abas 1 (Monitoramento), 2 (Tratamento), 3 (Modelo 3D), 4 (Cenários) e 5 (Engenharia): funcionais.
#
# Base dos resultados: simulação `titular` do modelo.jl (noradrenalina de fundo fixa), a mesma do relatório.
# O desmame da noradrenalina (`troca`) é opcional e usa premissas do grupo (aviso na tela).
# Nada do modelo.jl foi alterado: todos os números vêm de lá.
#
# Uso:  julia app_abas.jl
using GLMakie
using Printf

include(joinpath(@__DIR__, "modelo.jl"))
include(joinpath(@__DIR__, "aba_tratamento.jl"))
include(joinpath(@__DIR__, "aba_modelo3d.jl"))
include(joinpath(@__DIR__, "aba_engenharia.jl"))

# ================================================================ visual

const FUNDO    = RGBf(0.07, 0.09, 0.13)
const PAINEL   = RGBf(0.11, 0.14, 0.20)
const CARTAO   = RGBf(0.13, 0.17, 0.25)
const BORDA    = RGBf(0.22, 0.28, 0.38)
const TEXTO    = RGBf(0.90, 0.93, 0.97)
const GRIS     = RGBf(0.60, 0.65, 0.72)
const AZUL     = RGBf(0.38, 0.64, 1.00)    # noradrenalina
const LARANJA  = RGBf(1.00, 0.55, 0.30)    # angiotensina II
const VERDE    = RGBf(0.25, 0.82, 0.60)    # referência / sem BRA
const VERMELHO = RGBf(1.00, 0.35, 0.38)    # atenção / BRA
const AMARELO  = RGBf(0.97, 0.80, 0.30)    # IECA e "abaixo da meta"
const ATIVA    = AZUL
const INATIVA  = RGBf(0.20, 0.25, 0.34)

const USAR_SEGOE = true                     # troque para false se a fonte Segoe UI causar problema

const NOMES   = ["1 · Monitoramento", "2 · Tratamento", "3 · Modelo 3D", "4 · Cenários", "5 · Engenharia"]
const ROTULOS = ["Adulto, sem IECA/BRA", "Adulto, IECA", "Adulto, BRA",
                 "Idoso, sem IECA/BRA", "Idoso, IECA", "Idoso, BRA"]
const CURTOS  = ["A", "A-IECA", "A-BRA", "I", "I-IECA", "I-BRA"]
const USOS    = ["nenhum", "IECA", "BRA"]
# cor = uso crônico; traço = idade (tracejado: adulto, contínuo: idoso)
const COR     = [VERDE, AMARELO, VERMELHO, VERDE, AMARELO, VERMELHO]
const ESTILO  = [:dash, :dash, :dash, :solid, :solid, :solid]
const LARG    = [2.5, 5.0, 3.0, 2.5, 5.0, 3.0]

const R_VARREDURA   = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20]   # varredura de r (1 a 20)
const PASSO_GRAFICO = 10                                           # 1 ponto a cada 10 s

# premissas do grupo para a noradrenalina (mesmos valores iniciais dos sliders do app.jl)
const NE_FRAC = 0.5
const NE_EC50 = DOSE_NE_ALTA
const NE_T12  = T_MEIA_NE_FAIXA[2]

const AVISO = "Simulação para apoiar a discussão; não é diagnóstico. A farmacodinâmica da noradrenalina ainda não foi calibrada, " *
              "e a calibração de IECA/BRA vem de subgrupos pequenos (29 e 22 pacientes) de uma análise post hoc do ATHOS-3."
const AVISO_NE = "Premissas do grupo, não calibradas: Emax, EC50, t1/2 e regra de desmame da noradrenalina. Não usar como resultado clínico."

# ================================================================ utilidades

f0(x)  = @sprintf("%.0f", x)
f1(x)  = replace(@sprintf("%.1f", x), "." => ",")
f2(x)  = replace(@sprintf("%.2f", x), "." => ",")
fs1(x) = replace(@sprintf("%+.1f", x), "." => ",")
comalfa(c, a) = RGBAf(c.r, c.g, c.b, a)

function classe_f15(f)
    f > F15_CORTES[2] && return "responde"
    f < F15_CORTES[1] && return "não responde"
    return "zona cinzenta"
end

# cor do cartão da PAM: verde na meta, amarelo entre 65 e a meta, vermelho abaixo de 65
function cor_pam(p)
    p >= PAM_ALVO_ATHOS - 0.5 && return VERDE
    p >= PAM_ALVO_SSC && return AMARELO
    return VERMELHO
end

# ================================================================ estado compartilhado entre as abas

const PARAM = (Emax    = Observable(EMAX_CENTRAL),
               ganho   = Observable(1.0),
               fc      = Observable(FC_CHOQUE),
               dci     = Observable(round(DC_CHOQUE_IDOSO, digits = 1)),
               frp     = Observable(FATOR_RP_IDADE),
               fcx     = Observable(FATOR_C_IDADE),
               ne_frac = Observable(NE_FRAC),
               ne_ec50 = Observable(NE_EC50),
               ne_t12  = Observable(NE_T12))

const PAC     = Observable(2)                 # 1 adulto, 2 idoso
const USO     = Observable(3)                 # 1 nenhum, 2 IECA, 3 BRA
const K       = lift((p, u) -> 3 * (p - 1) + u, PAC, USO)    # índice do cenário (1 a 6)
const FONTE   = Observable(:titular)          # :titular (nora fixa) ou :troca (desmame da nora, premissas do grupo)
const SEL     = [Observable(k in (1, 4, 6)) for k in 1:6]    # cenários exibidos na aba 4

const DADOS   = Observable{Any}(nothing)      # resultados de todas as simulações
const INST    = Observable{Any}(nothing)      # estado no instante (minuto inteiro) atual
const TSIM    = Observable(0.0)               # min desde o início da Ang II
const TOCANDO = Observable(true)
const VEL     = Observable(3.0)               # min simulados por segundo
const IMIN    = Ref(0)
const PROG    = Ref(false)                    # true quando o programa (e não o usuário) move um slider
const ABA     = Ref(1)
const VERS    = Ref(0)                        # conta quantas vezes o modelo foi recalculado
const ATUALIZAR_SENS = Ref{Any}(() -> nothing)   # preenchida pela aba 4
const TBEAT   = Observable(0.0)               # relógio do batimento (s), só para animar a onda de pressão
const BEATS   = Dict{Tuple{Int,Float64},Any}()   # batimentos já calculados (chave: cenário e R_p)

# ================================================================ cálculo (só chama o modelo.jl)

function calcular(Emax, ganho, fc, dci, frp, fcx, ne_frac, ne_ec50, ne_t12)
    adulto  = planta_choque_adulto(FC = fc)
    idoso   = planta_choque_idoso(FC = fc, DC = dci, fRp = frp, fC = fcx)
    plantas = (adulto, adulto, adulto, idoso, idoso, idoso)
    rb      = r_bra(Emax)
    rs      = (1.0, 1.0, rb, 1.0, 1.0, rb)
    EC50    = ec50_aii(Emax)
    bas     = [pam_regime(plantas[k]) for k in 1:6]
    sims    = [titular(plantas[k], Emax, rs[k]; ganho) for k in 1:6]
    trocas  = [troca(plantas[k], Emax, rs[k]; ganho, frac_emax_ne = ne_frac, ec50_ne = ne_ec50, t_meia_ne = ne_t12)
               for k in 1:6]

    # séries no mesmo formato; tempo em min desde o início da Ang II
    series = Dict{Symbol,Any}(
        :titular => [(t = sims[k].t, pam = sims[k].pam, aii = sims[k].dose, nora = nothing,
                      x_aii = sims[k].x, x_ne = nothing,
                      e_aii = hill.(sims[k].x, Emax, rs[k] * EC50), e_ne = nothing) for k in 1:6],
        :troca   => [(t = trocas[k].t .- FASE_NE, pam = trocas[k].pam, aii = trocas[k].aii, nora = trocas[k].nora,
                      x_aii = trocas[k].aii, x_ne = trocas[k].nora,
                      e_aii = trocas[k].efeito_aii, e_ne = trocas[k].efeito_nora) for k in 1:6])

    # batimentos (Euler 0,1 ms): basal e ao fim das 3 h; IECA = sem uso (r = 1), então não repete a conta
    b0 = Vector{Any}(undef, 6)
    b0[1] = batimento(plantas[1]);  b0[2] = b0[1];  b0[3] = b0[1]
    b0[4] = batimento(plantas[4]);  b0[5] = b0[4];  b0[6] = b0[4]
    b3 = Dict{Symbol,Any}()
    for f in (:titular, :troca)
        v = Vector{Any}(undef, 6)
        for k in (1, 3, 4, 6)
            pam = series[f][k].pam[end]
            v[k] = batimento(com_rp(plantas[k], plantas[k].Rp * pam / bas[k]))
        end
        v[2] = v[1];  v[5] = v[4]
        b3[f] = v
    end

    f15 = [[titular(p, Emax, r; ganho) for r in R_VARREDURA] for p in (adulto, idoso)]
    return (; Emax, ganho, rb, rs, plantas, bas, sims, trocas, series, b0, b3, f15, EC50)
end

function recalcular!()
    VERS[] += 1
    empty!(BEATS)
    DADOS[] = calcular(PARAM.Emax[], PARAM.ganho[], PARAM.fc[], PARAM.dci[], PARAM.frp[], PARAM.fcx[],
                       PARAM.ne_frac[], PARAM.ne_ec50[], PARAM.ne_t12[])
    return
end

serie(d, k) = d.series[FONTE[]][k]
indice(s, tmin) = clamp(round(Int, (tmin - s.t[1]) * 60) + 1, 1, length(s.t))

# tempo até a meta, contado desde o início da Ang II (mesmo critério do modelo: PAM que arredonda para 75)
function t_alvo_serie(s)
    i0 = indice(s, 0.0)
    j  = findfirst(>=(PAM_ALVO_ATHOS - 0.5), @view s.pam[i0:end])
    return j === nothing ? nothing : s.t[i0 + j - 1]
end

# batimento (Euler 0,1 ms) para o R_p do instante; guarda o resultado para não refazer a conta quando
# o R_p não mudou (no platô da titulação, por exemplo). A chave arredonda R_p na 5ª casa decimal.
function batimento_cache(d, j, pam)
    Rp = d.plantas[j].Rp * pam / d.bas[j]
    return get!(BEATS, (j, round(Rp, digits = 5))) do
        batimento(com_rp(d.plantas[j], Rp))
    end
end

# estado no minuto inteiro atual: PAM, doses e ondas de pressão (só os batimentos que a tela usa)
function atualizar_instante!()
    d = DADOS[]
    d === nothing && return
    k  = K[]
    tm = Float64(IMIN[])
    bat  = Dict{Int,Any}()
    pams = Dict{Int,Float64}()
    for j in unique([k, 1, 4, 6])
        s   = serie(d, j)
        pam = s.pam[indice(s, tm)]
        pams[j] = pam
        bat[j]  = batimento_cache(d, j, pam)
    end
    s = serie(d, k);  i = indice(s, tm)
    INST[] = (; k, tm, bat, pam = pams[k], aii = s.aii[i], nora = s.nora === nothing ? nothing : s.nora[i],
              x_aii = s.x_aii[i], x_ne = s.x_ne === nothing ? nothing : s.x_ne[i],
              e_aii = s.e_aii[i], e_ne = s.e_ne === nothing ? nothing : s.e_ne[i])
    return
end

# ================================================================ textos das tabelas

# tabela da aba 4 (colunas compactadas e rodapé em 2 linhas para não cortar na margem direita)
function tabela_cenarios(d, f)
    cab = rpad("Cenário", 20) * lpad("PAM 0", 7) * lpad("PAS/PAD 0", 11) * lpad("PAM 3h", 8) *
          lpad("PAS/PAD 3h", 12) * lpad("Ang II", 8) * lpad("Nora", 9) * lpad("Meta 75", 10) *
          lpad("F15", 6) * "  " * "Classe F15 / status"
    linhas = String[]
    for k in 1:6
        SEL[k][] || continue
        s  = d.series[f][k]
        ta = t_alvo_serie(s)
        b0 = d.b0[k];  b3 = d.b3[f][k]
        nora = s.nora === nothing ? "fixa" : f2(s.nora[end]) * " µg"
        f15  = f == :titular ? f2(d.sims[k].F15) : "-"
        cls  = f == :titular ? classe_f15(d.sims[k].F15) : "-"
        status = ta === nothing ? "não atingida" : "atingida"
        push!(linhas, rpad(ROTULOS[k], 20) * lpad(f1(d.bas[k]), 7) * lpad(f0(b0.PAS) * "/" * f0(b0.PAD), 11) *
                      lpad(f1(s.pam[end]), 8) * lpad(f0(b3.PAS) * "/" * f0(b3.PAD), 12) *
                      lpad(f0(s.aii[end]), 8) * lpad(nora, 9) *
                      lpad(ta === nothing ? "não" : f0(ta) * " min", 10) * lpad(f15, 6) * "  " * cls * " / " * status)
    end
    isempty(linhas) && push!(linhas, "Nenhum cenário selecionado.")
    rod1 = "Tempos desde o início da Ang II (mmHg; Ang II em ng/kg/min). IECA = sem uso no modelo (r = 1)."
    rod2 = f == :titular ? "Noradrenalina de fundo fixa (não modelada)." :
                           "Com desmame da nora (premissas do grupo); F15 só existe na simulação base."
    return join(vcat([cab], linhas, ["", rod1, rod2]), "\n")
end

# sensibilidade: varia um parâmetro por vez (compactado para caber na coluna 4 sem cortar)
function texto_sens(par)
    d     = DADOS[]
    nomes = ("R_p com a idade (x)", "C com a idade (x)", "Emax da Ang II", "Fator r do BRA")
    vals  = (collect(1.10:0.02:1.20), collect(0.25:0.05:0.50), [0.3, 0.4, 0.5, 0.6, 0.8, 1.0],
             [1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0])[par]
    bra   = USO[] == 3
    l = String["Sensibilidade: " * nomes[par],
               "Idoso · uso crônico: " * USOS[USO[]] * " (nora fixa)",
               lpad("valor", 6) * lpad("PAM 3h", 8) * lpad("Ang II", 8) * lpad("Meta 75", 9) * lpad("PAS/PAD", 10)]
    for v in vals
        Emx = par == 3 ? v : PARAM.Emax[]
        r   = par == 4 ? v : (bra ? r_bra(Emx) : 1.0)
        frp = par == 1 ? v : PARAM.frp[]
        fcx = par == 2 ? v : PARAM.fcx[]
        p = planta_choque_idoso(FC = PARAM.fc[], DC = PARAM.dci[], fRp = frp, fC = fcx)
        s = titular(p, Emx, r; ganho = PARAM.ganho[])
        b = batimento(com_rp(p, s.rp_final))
        push!(l, lpad(f2(v), 6) * lpad(f1(s.pam_3h), 8) * lpad(f0(s.dose_final), 8) *
                 lpad(s.t_alvo === nothing ? "não" : f0(s.t_alvo) * "m", 9) *
                 lpad(f0(b.PAS) * "/" * f0(b.PAD), 10))
    end
    push!(l, "")
    par == 2 && push!(l, "C não entra na PAM média:\nmuda só PAS/PAD.")
    par == 3 && bra && push!(l, "r do BRA recalibrado\npara cada Emax.")
    par == 4 && push!(l, "r multiplica o EC50\n(EC50' = r x EC50).")
    push!(l, "Demais parâmetros nos\nvalores atuais.")
    return join(l, "\n")
end

# ================================================================ peças de interface

function barra!(gl, ativa, trocar!)
    for j in 1:5
        bt = Button(gl[1, j]; label = NOMES[j], width = 200,
                    buttoncolor = j == ativa ? ATIVA : INATIVA,
                    labelcolor = j == ativa ? :black : :white)
        on(_ -> trocar!(j), bt.clicks)
    end
    return
end

# botão escuro com texto claro (no tema escuro o padrão fica branco sobre cinza-claro e some)
function botao(pos, rotulo; kw...)
    try
        return Button(pos; label = rotulo, buttoncolor = INATIVA, labelcolor = :white,
                      buttoncolor_hover = BORDA, labelcolor_hover = :white, kw...)
    catch
        return Button(pos; label = rotulo, buttoncolor = INATIVA, labelcolor = :white, kw...)
    end
end

# controles de reprodução (iguais em todas as abas, ligados ao mesmo estado)
function controles!(gl)
    bp = botao(gl[1, 1], TOCANDO[] ? "Pausar" : "Reproduzir"; width = 120)
    bz = botao(gl[1, 2], "Voltar a 0"; width = 120)
    Label(gl[1, 3], "Instante"; halign = :right)
    st = Slider(gl[1, 4]; range = 0:1:Int(HORIZONTE), startvalue = floor(Int, TSIM[]), width = 420)
    Label(gl[1, 5], lift(t -> f0(floor(t)) * " min", TSIM); halign = :left, width = 70)
    Label(gl[1, 6], "Velocidade"; halign = :right)
    sv = Slider(gl[1, 7]; range = 0.5:0.5:10, startvalue = VEL[], width = 200)
    Label(gl[1, 8], lift(v -> f1(v) * " min/s", VEL); halign = :left, width = 90)

    on(bp.clicks) do _
        if !TOCANDO[] && TSIM[] >= HORIZONTE
            TSIM[] = 0.0
        end
        TOCANDO[] = !TOCANDO[]
    end
    on(TOCANDO) do p
        bp.label[] = p ? "Pausar" : "Reproduzir"
    end
    on(bz.clicks) do _
        TSIM[] = 0.0
    end
    on(TSIM) do t
        PROG[] = true
        set_close_to!(st, floor(t))
        PROG[] = false
    end
    on(st.value) do v
        PROG[] && return
        TSIM[] = Float64(v)
    end
    on(sv.value) do v
        v != VEL[] && (VEL[] = v)
    end
    on(VEL) do v
        sv.value[] != v && set_close_to!(sv, v)
    end
    return
end

# cartão de indicador: título pequeno, valor grande, linha de apoio
function cartao!(gl, titulo, valor, apoio, cor)
    try
        Box(gl[1:3, 1]; color = CARTAO, strokecolor = BORDA, strokewidth = 1.5, cornerradius = 12)
    catch
        Box(gl[1:3, 1]; color = CARTAO, strokecolor = BORDA, strokewidth = 1.5)    # versões sem cantos arredondados
    end
    Label(gl[1, 1], titulo; fontsize = 13, color = GRIS, tellwidth = false, padding = (0, 0, 0, 8))
    Label(gl[2, 1], valor; fontsize = 30, font = :bold, color = cor, tellwidth = false)
    Label(gl[3, 1], apoio; fontsize = 12, color = GRIS, tellwidth = false, padding = (0, 0, 8, 0))
    rowgap!(gl, 2)
    return
end

# ================================================================ aba 1: monitoramento hemodinâmico

function aba_monitoramento(corpo)
    # ---------------------------------------------------------------- perfil do paciente
    painel = corpo[1, 1] = GridLayout(tellheight = false, valign = :top)
    Label(painel[1, 1], "Perfil do paciente"; fontsize = 18, font = :bold, halign = :left)
    Label(painel[2, 1], "Idade"; halign = :left, color = GRIS)
    m_pac = Menu(painel[3, 1]; options = ["Adulto (46 anos)", "Idoso (70 anos)"], default = "Idoso (70 anos)")
    Label(painel[4, 1], "Uso crônico"; halign = :left, color = GRIS)
    m_uso = Menu(painel[5, 1]; options = ["Nenhum", "IECA", "BRA"], default = "BRA")
    on(m_pac.i_selected) do i
        i === nothing || (PAC[] = i)
    end
    on(m_uso.i_selected) do i
        i === nothing || (USO[] = i)
    end
    Label(painel[6, 1], "Parâmetros do choque"; halign = :left, color = GRIS)
    sg = SliderGrid(painel[7, 1],
        (label = "FC (bpm)", range = 63:1:120, startvalue = PARAM.fc[], format = "{:.0f}"),
        (label = "DC do idoso (L/min)", range = 7:0.1:10.2, startvalue = PARAM.dci[], format = "{:.1f}"),
        (label = "R_p com a idade (×)", range = 1.10:0.01:1.20, startvalue = PARAM.frp[], format = "{:.2f}"),
        (label = "C com a idade (×)", range = 0.25:0.01:0.50, startvalue = PARAM.fcx[], format = "{:.2f}"))
    for (s, o) in zip(sg.sliders, (PARAM.fc, PARAM.dci, PARAM.frp, PARAM.fcx))
        on(s.value) do v
            o[] != v && (o[] = v)
        end
        on(o) do v
            s.value[] != v && set_close_to!(s, v)
        end
    end
    Label(painel[8, 1], "Meta de PAM: " * f0(PAM_ALVO_ATHOS) * " mmHg (ATHOS-3)\nReferência: " * f0(PAM_ALVO_SSC) *
                        " mmHg (SSC)";
          halign = :left, justification = :left, fontsize = 13)
    lt = painel[9, 1] = GridLayout(halign = :left)
    tg = Toggle(lt[1, 1]; active = false)
    Label(lt[1, 2], "incluir desmame da noradrenalina\n(premissas do grupo)"; halign = :left, justification = :left)
    on(tg.active) do a
        FONTE[] = a ? :troca : :titular
    end
    on(FONTE) do f
        tg.active[] != (f == :troca) && (tg.active[] = (f == :troca))
    end
    Label(painel[10, 1], lift(f -> f == :troca ? AVISO_NE : "", FONTE);
          halign = :left, justification = :left, color = VERMELHO, fontsize = 12, word_wrap = true, width = 330)
    Label(painel[11, 1], AVISO; halign = :left, justification = :left, color = GRIS, fontsize = 11, word_wrap = true, width = 330)
    rowgap!(painel, 6)

    # ---------------------------------------------------------------- cartões
    direita = corpo[1, 2] = GridLayout()
    cards   = direita[1, 1] = GridLayout()
    defs = [
        ("PAM atual", lift(i -> f1(i.pam), INST),
         lift(i -> "mmHg · meta " * f0(PAM_ALVO_ATHOS) * " (" * fs1(i.pam - PAM_ALVO_ATHOS) * ")", INST),
         lift(i -> cor_pam(i.pam), INST)),
        ("PAS", lift(i -> f0(i.bat[i.k].PAS), INST), "mmHg · sistólica", TEXTO),
        ("PAD", lift(i -> f0(i.bat[i.k].PAD), INST), "mmHg · diastólica", TEXTO),
        ("Pressão de pulso", lift(i -> f0(i.bat[i.k].PP), INST), "mmHg · PAS menos PAD", TEXTO),
        ("Angiotensina II", lift(i -> f1(i.aii), INST), "ng/kg/min", LARANJA),
        ("Noradrenalina", lift(i -> i.nora === nothing ? "fixa" : f2(i.nora), INST),
         lift(i -> i.nora === nothing ? "não modelada" : "µg/kg/min", INST), AZUL)]
    for (j, (tit, val, apoio, cor)) in enumerate(defs)
        sub = cards[1, j] = GridLayout()
        cartao!(sub, tit, val, apoio, cor)
    end
    colgap!(cards, 10)

    # ---------------------------------------------------------------- gráficos
    gg = direita[2, 1] = GridLayout()

    # 1. evolução da PAM
    ax_pam = Axis(gg[1, 1]; title = "Evolução da PAM", xlabel = "tempo desde o início da Ang II (min)",
                  ylabel = "PAM (mmHg)", xticks = -60:30:180)
    hlines!(ax_pam, [PAM_ALVO_ATHOS]; color = VERDE, linestyle = :dash, linewidth = 1.5)
    hlines!(ax_pam, [PAM_ALVO_SSC]; color = GRIS, linestyle = :dash, linewidth = 1.5)
    text!(ax_pam, HORIZONTE, PAM_ALVO_ATHOS; text = "meta (75) ", align = (:right, :bottom), color = VERDE, fontsize = 12)
    text!(ax_pam, HORIZONTE, PAM_ALVO_SSC; text = "65 ", align = (:right, :top), color = GRIS, fontsize = 12)
    obs_pam = Observable(Point2f[(0, 0)])
    cor_sel = lift(k -> COR[k], K)
    lines!(ax_pam, obs_pam; color = cor_sel, linewidth = 3.5)
    vlines!(ax_pam, TSIM; color = (:white, 0.45))
    marca = Observable([Point2f(0, 0)])
    scatter!(ax_pam, marca; color = :white, markersize = 11)

    # 2. forma de onda da pressão arterial
    sub_onda = Observable("")
    ax_onda = Axis(gg[1, 2]; title = "Forma de onda da pressão arterial", subtitle = sub_onda, subtitlesize = 12,
                   ylabel = "pressão (mmHg)", xlabel = "tempo (s)")
    obs_onda = [Observable(Point2f[(0, 0)]) for _ in 1:3]
    for (j, kk) in enumerate((1, 4, 6))
        lines!(ax_onda, obs_onda[j]; color = COR[kk], linestyle = ESTILO[kk], linewidth = 2, label = ROTULOS[kk])
    end
    onda_sel = Observable(Point2f[(0, 0)])
    lines!(ax_onda, onda_sel; color = :white, linewidth = 1.5, linestyle = :dot, label = "cenário selecionado")
    pas_pad = Observable([0.0, 0.0])
    hlines!(ax_onda, pas_pad; color = (:white, 0.35), linestyle = :dash, linewidth = 1)
    axislegend(ax_onda; position = :rt, framevisible = false, labelsize = 11, patchsize = (16, 10), rowgap = 0)
    marca_onda = Observable([Point2f(0, 0)])
    scatter!(ax_onda, marca_onda; color = VERMELHO, markersize = 12)

    # 3. vasopressores: Ang II (eixo esquerdo) e nora (eixo direito, só no desmame)
    sub_d = Observable("")
    ax_d = Axis(gg[1, 3]; title = "Administração dos vasopressores", subtitle = sub_d, subtitlesize = 12,
                xlabel = "tempo desde o início da Ang II (min)", xticks = -60:30:180,
                ylabel = "Ang II (ng/kg/min)", ylabelcolor = LARANJA, yticklabelcolor = LARANJA)
    ax_n = Axis(gg[1, 3]; yaxisposition = :right, ylabel = "nora (µg/kg/min)", ylabelcolor = AZUL,
                yticklabelcolor = AZUL, backgroundcolor = :transparent)
    hidexdecorations!(ax_n);  hidespines!(ax_n);  hideydecorations!(ax_n; label = false, ticklabels = false, ticks = false)
    linkxaxes!(ax_d, ax_n)
    hlines!(ax_d, [DOSE_MAX]; color = GRIS, linestyle = :dash, linewidth = 1)
    obs_aii  = Observable(Point2f[(0, 0)])
    obs_nora = Observable(Point2f[(NaN, NaN)])
    lines!(ax_d, obs_aii; color = LARANJA, linewidth = 3)
    lines!(ax_n, obs_nora; color = AZUL, linewidth = 3)
    vlines!(ax_d, TSIM; color = (:white, 0.45))
    marca_aii = Observable([Point2f(0, 0)])
    scatter!(ax_d, marca_aii; color = LARANJA, markersize = 11)

    function desenhar_curvas!()
        d = DADOS[];  f = FONTE[];  k = K[]
        s = d.series[f][k]
        ii = 1:PASSO_GRAFICO:length(s.t)
        x0 = f == :troca ? -FASE_NE : 0.0
        obs_pam[] = Point2f.(s.t[ii], s.pam[ii])
        obs_aii[] = Point2f.(s.t[ii], s.aii[ii])
        obs_nora[] = s.nora === nothing ? [Point2f(NaN, NaN)] : Point2f.(s.t[ii], s.nora[ii])
        limits!(ax_pam, x0, HORIZONTE, min(55.0, floor(minimum(s.pam) - 1)), max(80.0, ceil(maximum(s.pam) + 1)))
        limits!(ax_d, x0, HORIZONTE, 0, 85)
        limits!(ax_n, x0, HORIZONTE, 0, DOSE_NE_MAX * 1.1)
        vis = f == :troca                      # o eixo da nora só aparece quando há curva de nora
        ax_n.ylabelvisible[] = vis
        ax_n.yticklabelsvisible[] = vis
        ax_n.yticksvisible[] = vis
        sub_d[] = f == :troca ? "nora: premissas do grupo (azul, eixo direito)" : "noradrenalina: dose de fundo fixa"
        return
    end

    function desenhar_marcas!()
        d = DADOS[];  s = serie(d, K[]);  i = indice(s, TSIM[])
        marca[]     = [Point2f(TSIM[], s.pam[i])]
        marca_aii[] = [Point2f(TSIM[], s.aii[i])]
        return
    end

    function desenhar_ondas!()
        inst = INST[];  d = DADOS[]
        ymax = 0.0
        for (j, kk) in enumerate((1, 4, 6))
            b = inst.bat[kk];  T = 60 / d.plantas[kk].FC
            obs_onda[j][] = vcat(Point2f.(b.t, b.pao), Point2f.(b.t .+ T, b.pao), Point2f.(b.t .+ 2T, b.pao))
            ymax = max(ymax, b.PAS)
        end
        b = inst.bat[inst.k];  T = 60 / d.plantas[inst.k].FC
        onda_sel[] = vcat(Point2f.(b.t, b.pao), Point2f.(b.t .+ T, b.pao), Point2f.(b.t .+ 2T, b.pao))
        pas_pad[]  = [b.PAD, b.PAS]
        sub_onda[] = ROTULOS[inst.k] * ": " * f0(b.PAS) * "/" * f0(b.PAD) * " mmHg  (t = " * f0(inst.tm) * " min)"
        limits!(ax_onda, 0, 3 * T, 0, ymax * 1.5)
        return
    end

    on(_ -> (desenhar_curvas!(); desenhar_marcas!()), DADOS)
    on(_ -> (desenhar_curvas!(); desenhar_marcas!()), FONTE)
    on(_ -> (desenhar_curvas!(); desenhar_marcas!()), K)
    on(_ -> desenhar_marcas!(), TSIM)
    # ponto que percorre a onda de pressão, como num monitor (3 ciclos, em loop)
    function desenhar_marca_onda!()
        inst = INST[];  d = DADOS[]
        b = inst.bat[inst.k];  T = 60 / d.plantas[inst.k].FC
        x = mod(TBEAT[], 3T)
        i = clamp(round(Int, mod(x, T) / DT_EULER) + 1, 1, length(b.pao))
        marca_onda[] = [Point2f(x, b.pao[i])]
        return
    end

    on(_ -> desenhar_ondas!(), INST)
    on(_ -> desenhar_marca_onda!(), TBEAT)
    desenhar_curvas!();  desenhar_marcas!();  desenhar_ondas!();  desenhar_marca_onda!()
    colsize!(corpo, 1, Fixed(350))
    return
end

# ================================================================ aba 4: análise de cenários

function aba_cenarios(corpo)
    # ---------------------------------------------------------------- painel esquerdo
    painel = corpo[1, 1] = GridLayout(tellheight = false, valign = :top)
    Label(painel[1, 1], "Cenários exibidos"; fontsize = 18, font = :bold, halign = :left)
    for k in 1:6
        lin = painel[k + 1, 1] = GridLayout(halign = :left)
        tg = Toggle(lin[1, 1]; active = SEL[k][])
        Label(lin[1, 2], ROTULOS[k]; color = COR[k], halign = :left)
        on(tg.active) do a
            SEL[k][] = a
        end
    end
    Label(painel[8, 1], "Cor: uso crônico (verde nenhum, amarelo IECA, vermelho BRA).\nTracejado: adulto. Contínuo: idoso.";
          halign = :left, justification = :left, color = GRIS, fontsize = 12)
    Label(painel[9, 1], "Parâmetros de entrada"; halign = :left, font = :bold)
    sg = SliderGrid(painel[10, 1],
        (label = "Emax da Ang II (0,3 a 1,0)", range = 0.3:0.05:1.0, startvalue = PARAM.Emax[], format = "{:.2f}"),
        (label = "Passo: ng/kg/min por mmHg", range = 0.3:0.1:2.0, startvalue = PARAM.ganho[], format = "{:.1f}"))
    for (s, o) in zip(sg.sliders, (PARAM.Emax, PARAM.ganho))
        on(s.value) do v
            o[] != v && (o[] = v)
        end
        on(o) do v
            s.value[] != v && set_close_to!(s, v)
        end
    end
    Label(painel[11, 1], lift(d -> "Calibração (Leisman 2024): EC50 = " * f1(d.EC50) * " ng/kg/min\nr do BRA para este Emax = " * f1(d.rb), DADOS);
          halign = :left, justification = :left, color = GRIS, fontsize = 12)
    Label(painel[12, 1], "Análise de sensibilidade (parâmetro variado)"; halign = :left, font = :bold)
    m_sens = Menu(painel[13, 1]; options = ["R_p com a idade", "C com a idade", "Emax da Ang II", "Fator r do BRA"],
                  default = "Fator r do BRA")
    Label(painel[14, 1], "EC50 da Ang II entra por r (EC50' = r x EC50).\nEC50 da noradrenalina: aba 2 (premissa).";
          halign = :left, justification = :left, color = GRIS, fontsize = 12)
    Label(painel[15, 1], lift(f -> f == :troca ? AVISO_NE : "", FONTE);
          halign = :left, justification = :left, color = VERMELHO, fontsize = 12, word_wrap = true, width = 330)
    Label(painel[16, 1], "Paciente e uso crônico: escolha na aba 1.\nDesmame da nora: ligue na aba 1 ou 2.";
          halign = :left, justification = :left, color = GRIS, fontsize = 12)
    rowgap!(painel, 6)

    # ---------------------------------------------------------------- tabela e gráficos
    g = corpo[1, 2] = GridLayout()
    tab = Observable("")
    Label(g[1, 1:4], tab; font = "Consolas", fontsize = 12, halign = :left, justification = :left, tellwidth = false)

    # linha A: PAM, dose, F15, PAS/PAD
    ax_pam = Axis(g[2, 1]; title = "PAM", ylabel = "PAM (mmHg)", xlabel = "min desde o início da Ang II", xticks = -60:30:180)
    ax_dose = Axis(g[2, 2]; title = "Dose de Ang II", ylabel = "ng/kg/min", xlabel = "min desde o início da Ang II",
                   xticks = -60:30:180)
    hlines!(ax_pam, [PAM_ALVO_ATHOS]; color = VERDE, linestyle = :dash, linewidth = 1.5)
    hlines!(ax_pam, [PAM_ALVO_SSC]; color = GRIS, linestyle = :dash, linewidth = 1.5)
    hlines!(ax_dose, [DOSE_MAX]; color = GRIS, linestyle = :dash, linewidth = 1.5)
    obs_pam  = [Observable(Point2f[(0, 0)]) for _ in 1:6]
    obs_dose = [Observable(Point2f[(0, 0)]) for _ in 1:6]
    for k in (2, 5, 1, 4, 3, 6)     # IECA primeiro (mais grosso, por baixo): coincide com "nenhum"
        lines!(ax_pam, obs_pam[k]; color = COR[k], linestyle = ESTILO[k], linewidth = LARG[k], visible = SEL[k])
        lines!(ax_dose, obs_dose[k]; color = COR[k], linestyle = ESTILO[k], linewidth = LARG[k], visible = SEL[k])
    end
    vlines!(ax_pam, TSIM; color = (:white, 0.45));  vlines!(ax_dose, TSIM; color = (:white, 0.45))

    ax_f15 = Axis(g[2, 3]; title = "F15 (déficit fechado em 15 min)", subtitle = "círculo cheio: atinge 75 em 3 h",
                  subtitlesize = 11, xlabel = "fator r (log)", ylabel = "F15", xscale = log10,
                  xticks = ([1, 2, 3, 5, 10, 20], ["1", "2", "3", "5", "10", "20"]), limits = ((0.9, 22), (0, 1)))
    hspan!(ax_f15, F15_CORTES...; color = (:gray50, 0.3))
    text!(ax_f15, 21, F15_CORTES[2]; text = "zona cinzenta", align = (:right, :bottom), color = GRIS, fontsize = 12)
    obs_f15 = [Observable(Point2f[(1, 0)]) for _ in 1:2]
    cor_f15 = [Observable(fill(comalfa(PAINEL, 1), length(R_VARREDURA))) for _ in 1:2]
    cores_f15 = (RGBf(0.75, 0.80, 0.90), VERDE)
    for (j, nome) in enumerate(("Adulto", "Idoso"))
        lines!(ax_f15, obs_f15[j]; color = cores_f15[j], linewidth = 2.5, linestyle = j == 1 ? :dash : :solid)
        scatter!(ax_f15, obs_f15[j]; color = cor_f15[j], strokecolor = cores_f15[j], strokewidth = 2.5, markersize = 12)
        text!(ax_f15, 0.97, 0.97 - 0.07 * (j - 1); text = nome, space = :relative, align = (:right, :top),
              color = cores_f15[j], fontsize = 12)
    end
    r_atual = Observable(1.0)
    vlines!(ax_f15, r_atual; color = VERMELHO, linestyle = :dash)

    ax_pp = Axis(g[2, 4]; title = "PAS e PAD em 3 h", subtitle = "barra: PAD a PAS; losango: PAM", subtitlesize = 11,
                 ylabel = "mmHg", xticks = (1:6, CURTOS), xticklabelrotation = pi / 4)
    obs_pp = Observable(Point2f.(1:6, zeros(6)));  obs_pad = Observable(zeros(6))
    obs_ppcor = Observable([comalfa(COR[k], 0) for k in 1:6])
    barplot!(ax_pp, obs_pp; fillto = obs_pad, color = obs_ppcor, width = 0.55)
    obs_dots = Observable(Point2f.(1:6, zeros(6)))
    obs_dotcor = Observable([RGBAf(1, 1, 1, 0) for k in 1:6])
    scatter!(ax_pp, obs_dots; color = obs_dotcor, marker = :diamond, markersize = 11)

    # linha B: resistência, complacência, superfície 3D, sensibilidade
    ax_rp = Axis(g[3, 1]; title = "Resistência periférica R_p", subtitle = "cinza: choque; cor: após 3 h", subtitlesize = 11,
                 ylabel = "mmHg·s/ml", xticks = (1:6, CURTOS), xticklabelrotation = pi / 4)
    obs_rp0 = Observable(Point2f.(1:6, zeros(6)));  obs_rp3 = Observable(Point2f.(1:6, zeros(6)))
    obs_rp0cor = Observable([comalfa(GRIS, 0) for k in 1:6]);  obs_rp3cor = Observable([comalfa(COR[k], 0) for k in 1:6])
    barplot!(ax_rp, obs_rp0; color = obs_rp0cor, width = 0.36)
    barplot!(ax_rp, obs_rp3; color = obs_rp3cor, width = 0.36)
    ax_c = Axis(g[3, 2]; title = "Complacência arterial C", subtitle = "constante (não muda com fármaco)",
                subtitlesize = 11, ylabel = "ml/mmHg", xticks = (1:6, CURTOS), xticklabelrotation = pi / 4)
    obs_c = Observable(Point2f.(1:6, zeros(6)));  obs_ccor = Observable([comalfa(COR[k], 0) for k in 1:6])
    barplot!(ax_c, obs_c; color = obs_ccor, width = 0.55)

    ax_3d = Axis3(g[3, 3]; title = "PAM em regime do idoso: dose × r\nplano: alvo 75; pontos: fim das 3 h (base)",
                  titlesize = 13, protrusions = (70, 30, 40, 45),
                  xlabel = "dose (ng/kg/min)", ylabel = "fator r", zlabel = "PAM (mmHg)",
                  yticks = (log10.([1, 2, 5, 10, 20]), ["1", "2", "5", "10", "20"]),
                  azimuth = 1.3pi, elevation = 0.12pi)
    doses = range(0, DOSE_MAX, length = 41)
    logr  = range(0, log10(20), length = 41)
    sup   = Observable(zeros(Float32, 41, 41))
    surface!(ax_3d, doses, logr, sup; colormap = :viridis, colorrange = (60, 85))
    surface!(ax_3d, [0, DOSE_MAX], [0, log10(20)], fill(Float32(PAM_ALVO_ATHOS), 2, 2);
             colormap = [(:white, 0.14), (:white, 0.14)], transparency = true)
    pts3d = Observable([Point3f(0), Point3f(0)])
    scatter!(ax_3d, pts3d; color = [VERDE, VERMELHO], markersize = 16, strokecolor = :white, strokewidth = 1.5)

    sens = Observable("")
    Label(g[3, 4], sens; font = "Consolas", fontsize = 11, halign = :left, valign = :top, justification = :left,
          tellwidth = false, tellheight = false)

    # ---------------------------------------------------------------- atualização
    function atualizar!()
        d = DADOS[];  f = FONTE[];  ser = d.series[f]
        x0   = f == :troca ? -FASE_NE : 0.0
        ymin = min(55.0, floor(minimum(minimum(s.pam) for s in ser) - 1))
        ymax = max(80.0, ceil(maximum(maximum(s.pam) for s in ser) + 1))
        for k in 1:6
            s = ser[k];  ii = 1:PASSO_GRAFICO:length(s.t)
            obs_pam[k][]  = Point2f.(s.t[ii], s.pam[ii])
            obs_dose[k][] = Point2f.(s.t[ii], s.aii[ii])
        end
        limits!(ax_pam, x0, HORIZONTE, ymin, ymax)
        limits!(ax_dose, x0, HORIZONTE, 0, 85)

        for (j, v) in enumerate(d.f15)
            cor = cores_f15[j]
            cor_f15[j][] = [s.t_alvo === nothing ? comalfa(PAINEL, 1) : comalfa(cor, 1) for s in v]
            obs_f15[j][] = Point2f.(R_VARREDURA, [s.F15 for s in v])
        end
        r_atual[] = d.rb

        sel = [SEL[k][] for k in 1:6]
        b3  = d.b3[f]
        obs_pp[]    = Point2f.(1:6, [b3[k].PAS for k in 1:6])
        obs_pad[]   = [b3[k].PAD for k in 1:6]
        obs_ppcor[] = [comalfa(COR[k], sel[k] ? 0.9 : 0.0) for k in 1:6]
        obs_dots[]  = Point2f.(1:6, [ser[k].pam[end] for k in 1:6])
        obs_dotcor[] = [RGBAf(1, 1, 1, sel[k] ? 1 : 0) for k in 1:6]
        limits!(ax_pp, 0.4, 6.6, 0, max(120.0, ceil(maximum(b3[k].PAS for k in 1:6) / 10) * 10 + 10))

        rp0 = [d.plantas[k].Rp for k in 1:6]
        rp3 = [d.plantas[k].Rp * ser[k].pam[end] / d.bas[k] for k in 1:6]
        obs_rp0[]    = Point2f.((1:6) .- 0.2, rp0)
        obs_rp3[]    = Point2f.((1:6) .+ 0.2, rp3)
        obs_rp0cor[] = [comalfa(GRIS, sel[k] ? 0.8 : 0.0) for k in 1:6]
        obs_rp3cor[] = [comalfa(COR[k], sel[k] ? 0.95 : 0.0) for k in 1:6]
        limits!(ax_rp, 0.4, 6.6, 0, maximum(vcat(rp0, rp3)) * 1.15)
        obs_c[]    = Point2f.(1:6, [d.plantas[k].C for k in 1:6])
        obs_ccor[] = [comalfa(COR[k], sel[k] ? 0.9 : 0.0) for k in 1:6]
        limits!(ax_c, 0.4, 6.6, 0, maximum(d.plantas[k].C for k in 1:6) * 1.15)

        bi = pam_regime(d.plantas[4])
        sup[] = Float32[pam_com_aii(bi, dd, d.Emax, d.EC50, 10.0^lr) for dd in doses, lr in logr]
        pts3d[] = [Point3f(d.sims[4].dose_final, 0, d.sims[4].pam_3h),
                   Point3f(d.sims[6].dose_final, log10(d.rb), d.sims[6].pam_3h)]
        ax_3d.limits[] = (0, DOSE_MAX, 0, log10(20), min(60.0, bi - 1), max(85.0, maximum(sup[]) + 1))

        tab[] = tabela_cenarios(d, f)
        return
    end

    # a sensibilidade refaz muitas simulações: só recalcula se algo mudou e se esta aba está à vista
    chave = Ref((0, 0, 0))
    function atualizar_sens!()
        c = (VERS[], something(m_sens.i_selected[], 4), USO[])
        c == chave[] && return
        chave[] = c
        sens[] = texto_sens(c[2])
        return
    end
    ATUALIZAR_SENS[] = atualizar_sens!

    on(_ -> (atualizar!(); ABA[] == 4 && atualizar_sens!()), DADOS)
    on(_ -> atualizar!(), FONTE)
    for k in 1:6
        on(_ -> atualizar!(), SEL[k])
    end
    on(_ -> atualizar_sens!(), m_sens.i_selected)
    atualizar!();  atualizar_sens!()
    colsize!(g, 1, Relative(0.23))
    colsize!(g, 2, Relative(0.23))
    colsize!(g, 3, Relative(0.24))
    colsize!(g, 4, Relative(0.30))
    colgap!(g, 14)
    colsize!(corpo, 1, Fixed(350))
    return
end

# ================================================================ tema e janela

function fonte_ou(nome, padrao)
    try
        GLMakie.Makie.to_font(nome)
        return nome
    catch
        return padrao
    end
end

function tema!()
    set_theme!(theme_dark())
    update_theme!(backgroundcolor = FUNDO, textcolor = TEXTO, Axis = (; backgroundcolor = PAINEL))
    if USAR_SEGOE
        try
            update_theme!(fonts = (; regular = fonte_ou("Segoe UI", "TeX Gyre Heros Makie"),
                                     bold = fonte_ou("Segoe UI Bold", "TeX Gyre Heros Makie Bold")))
        catch e
            @warn "Não consegui aplicar a fonte Segoe UI; usando a padrão." exception = e
        end
    end
    return
end

# avança o relógio da simulação (compartilhado por todas as abas)
function avancar!(dt)
    TBEAT[] = TBEAT[] + 0.12 * dt                # batimento a 0,4 do tempo real
    if TOCANDO[]
        t = TSIM[] + VEL[] * dt
        TSIM[] = t >= HORIZONTE ? t - HORIZONTE : t     # em loop: ao chegar em 180 min volta ao início
    end
    return
end

function main()
    tema!()
    recalcular!()
    atualizar_instante!()

    # o estado no minuto inteiro acompanha o relógio, o banco de dados e a escolha de cenário/fonte
    on(TSIM) do t
        m = floor(Int, t)
        if m != IMIN[]
            IMIN[] = m
            atualizar_instante!()
        end
    end
    on(_ -> atualizar_instante!(), DADOS)
    on(_ -> atualizar_instante!(), K)
    on(_ -> atualizar_instante!(), FONTE)
    onany((args...) -> recalcular!(), PARAM.Emax, PARAM.ganho, PARAM.fc, PARAM.dci, PARAM.frp, PARAM.fcx,
          PARAM.ne_frac, PARAM.ne_ec50, PARAM.ne_t12)

    figs = Vector{Any}(undef, 5)
    tela = Ref{Any}(nothing)
    function trocar!(j)
        ABA[] = j
        j == 4 && ATUALIZAR_SENS[]()
        display(tela[], figs[j])
        return
    end

    for i in 1:5
        fig = Figure(size = (1750, 980), backgroundcolor = FUNDO)
        gb = fig[1, 1] = GridLayout()
        barra!(gb, i, trocar!)
        corpo = fig[2, 1] = GridLayout()
        if i == 1
            aba_monitoramento(corpo)
        elseif i == 2
            aba_tratamento(corpo)
        elseif i == 3
            aba_modelo3d(corpo)
        elseif i == 4
            aba_cenarios(corpo)
        else
            aba_engenharia(corpo)
        end
        if i in (1, 2, 3, 4)
            ctrl = fig[3, 1] = GridLayout()
            controles!(ctrl)
        else
            rodape5 = fig[3, 1] = GridLayout()
            Label(rodape5[1, 1], "Documentação técnica e rastreabilidade sincronizadas em tempo real com modelo.jl";
                  fontsize = 12, color = GRIS, halign = :center)
        end
        rowsize!(fig.layout, 1, Fixed(44))
        rowsize!(fig.layout, 2, Relative(0.88))
        rowsize!(fig.layout, 3, Fixed(44))
        on(events(fig).tick) do tick
            avancar!(clamp(tick.delta_time, 0.0, 0.1))
        end
        figs[i] = fig
    end

    tela[] = display(figs[1])
    isinteractive() || wait(tela[])
    return
end

main()