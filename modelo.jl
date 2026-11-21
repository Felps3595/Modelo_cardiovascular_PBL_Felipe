# Choque vasodilatador refratário: Windkessel 0D de 3 elementos + angiotensina II (PK/PD)
# Todos os números vêm do documento "Choque vasodilatador refratário: Windkessel 0D com
# controle de vasopressores". Cada constante indica a seção/fonte de onde foi tirada.
#
# Unidades: pressão mmHg, fluxo ml/s, volume ml, tempo s (min onde indicado),
#           resistência mmHg·s/ml, complacência ml/mmHg, dose de Ang II ng/kg/min.

# ------------------------------------------------------------------ dados do documento

const P_VEN        = 0.0      # convenção P_ven = 0 (Melo & Peters, 1999)
const RAZAO_ZC     = 0.06     # Z_c / R_p (Westerhof 2009)
const DT_EULER     = 1e-4     # passo de Euler explícito: 0,1 ms

# Adulto saudável (46 anos): tabela "Entradas, parâmetros e saídas"
const FC_SAUDAVEL  = 63.0     # bpm (Alhakak 2023)
const RP_ADULTO    = 0.91     # ajuste para PA 124/74
const C_ADULTO     = 1.91     # ajuste para PA 124/74
vs_idade(idade)    = 106.4 - 0.27 * idade      # ml (Carlsson 2012)
const IDADE_ADULTO = 46.0
const IDADE_IDOSO  = 70.0

# Envelhecimento (Westerhof 2009): valores centrais e faixas do documento
const FATOR_RP_IDADE = 1.15   # faixa 1,10 a 1,20
const FATOR_C_IDADE  = 0.5    # faixa 0,25 a 0,50

# Choque (Melo & Peters 1999, sepse definitiva)
const RVS_CHOQUE   = 455.0    # dyn·s·cm⁻⁵, calculada como 80 × PAM / DC (sem pressão venosa)
const DC_CHOQUE_ADULTO = 10.19                       # L/min
const DC_CHOQUE_IDOSO  = 10.19 * (5.5 / 5.9)         # L/min = 9,5 (premissa do documento)
const FC_CHOQUE    = 90.0     # bpm (Bone 1992; premissa)

# Angiotensina II (bula GIAPREZA; protocolo ATHOS-3)
const T_MEIA_AII   = 1.0      # min (limite superior da bula)
const DOSE_INICIAL = 20.0     # ng/kg/min
const DOSE_MAX     = 80.0     # ng/kg/min nas 3 primeiras horas
const PASSO_MAX    = 15.0     # ng/kg/min por ajuste
const INTERVALO_AJUSTE = 5.0  # min
const PAM_ALVO_ATHOS   = 75.0 # mmHg
const PAM_ALVO_SSC     = 65.0 # mmHg
const HORIZONTE    = 180.0    # min

# Calibração da PD da Ang II (Leisman et al., 2024, subgrupos do ATHOS-3)
const CAL_SEM = (pam0 = 66.0, dose = 20.7, efeito = 9.1)          # efeito ajustado ao placebo
const CAL_BRA = (pam0 = 65.0, dose = 68.0, efeito = 5.3 - 2.2)    # +5,3 menos placebo +2,2
const EMAX_CENTRAL = 0.5
const EMAX_FAIXA   = (0.3, 1.0)

# Critério de não respondedor precoce
const F15_CORTES = (0.58, 0.61)

# ------------------------------------------------------------------ planta WK3

struct Planta
    nome::String
    FC::Float64     # bpm
    VS::Float64     # ml
    Rp::Float64     # mmHg·s/ml
    C::Float64      # ml/mmHg
end

t_ejecao(FC) = 0.413 - 0.0017 * FC                  # s (Weissler 1968)
zc(Rp)       = RAZAO_ZC * Rp
dc_mls(p::Planta)  = p.VS * p.FC / 60               # ml/s
dc_lmin(p::Planta) = p.VS * p.FC / 1000             # L/min

# RVS em dyn·s·cm⁻⁵ -> resistência total (R_p + Z_c) em mmHg·s/ml
rtotal_de_rvs(rvs) = rvs / 80 * 60 / 1000
const RP_CHOQUE_ADULTO = rtotal_de_rvs(RVS_CHOQUE) / (1 + RAZAO_ZC)     # 0,32

planta_adulto() = Planta("Adulto (46 anos)", FC_SAUDAVEL, vs_idade(IDADE_ADULTO), RP_ADULTO, C_ADULTO)
planta_idoso(; fRp = FATOR_RP_IDADE, fC = FATOR_C_IDADE) =
    Planta("Idoso (70 anos)", FC_SAUDAVEL, vs_idade(IDADE_IDOSO), RP_ADULTO * fRp, C_ADULTO * fC)
planta_choque_adulto(; FC = FC_CHOQUE) =
    Planta("Choque adulto", FC, DC_CHOQUE_ADULTO * 1000 / FC, RP_CHOQUE_ADULTO, C_ADULTO)
planta_choque_idoso(; FC = FC_CHOQUE, DC = DC_CHOQUE_IDOSO, fRp = FATOR_RP_IDADE, fC = FATOR_C_IDADE) =
    Planta("Choque idoso", FC, DC * 1000 / FC, RP_CHOQUE_ADULTO * fRp, C_ADULTO * fC)

com_rp(p::Planta, Rp) = Planta(p.nome, p.FC, p.VS, Rp, p.C)

# Fluxo de ejeção: meia-senoide com área VS e duração T_ej; zero no resto do ciclo
function q_in(t_no_ciclo, VS, Tej)
    return t_no_ciclo < Tej ? pi * VS / (2 * Tej) * sin(pi * t_no_ciclo / Tej) : 0.0
end

# PAM em regime periódico: média de P_ao = (R_p + Z_c) × DC + P_ven
pam_regime(p::Planta) = (1 + RAZAO_ZC) * p.Rp * dc_mls(p) + P_VEN

# Integra C dPc/dt = Q_in − (Pc − P_ven)/R_p por Euler explícito até o regime periódico
# e devolve o último batimento: P_ao(t) = Pc + Z_c Q_in
function batimento(p::Planta; dt = DT_EULER)
    T   = 60 / p.FC
    Tej = t_ejecao(p.FC)
    Zc  = zc(p.Rp)
    n   = round(Int, T / dt)
    nbat = ceil(Int, 8 * p.Rp * p.C / T) + 3
    Pc  = p.Rp * dc_mls(p) + P_VEN
    ts  = Vector{Float64}(undef, n)
    pao = Vector{Float64}(undef, n)
    for b in 1:nbat, i in 1:n
        tc = (i - 1) * dt
        q  = q_in(tc, p.VS, Tej)
        if b == nbat
            ts[i] = tc;  pao[i] = Pc + Zc * q
        end
        Pc += dt * (q - (Pc - P_VEN) / p.Rp) / p.C
    end
    pas, pad = maximum(pao), minimum(pao)
    return (t = ts, pao = pao, PAS = pas, PAD = pad, PAM = sum(pao) / n, PP = pas - pad)
end

# ------------------------------------------------------------------ angiotensina II

# PD de Hill com n = 1: E = Emax × x / (r × EC50 + x)
hill(x, Emax, EC50) = Emax * x / (EC50 + x)

# EC50 (em dose de regime, ng/kg/min) que reproduz o efeito ajustado sem IECA/BRA:
# +9,1 mmHg partindo de PAM 66 com 20,7 ng/kg/min
ec50_aii(Emax) = CAL_SEM.dose * (Emax * CAL_SEM.pam0 / CAL_SEM.efeito - 1)

# Fator r do BRA (Schild) que reproduz +3,1 mmHg partindo de PAM 65 com 68 ng/kg/min
r_bra(Emax) = (Emax * CAL_BRA.dose * CAL_BRA.pam0 / CAL_BRA.efeito - CAL_BRA.dose) / ec50_aii(Emax)

# PAM quase-estática: a constante de tempo do Windkessel (R_p × C ≈ 0,6 s no choque) é muito
# menor que a da PK (minutos), então PAM ≈ 1,06 × DC × R_p,choque × [1 + E_AII]
pam_com_aii(pam_basal, x, Emax, EC50, r) = pam_basal * (1 + hill(x, Emax, r * EC50))

# Maior PAM alcançável com a Ang II saturada (noradrenalina de fundo fixa, já na planta de choque)
pam_max(pam_basal, Emax, EC50, r) = pam_com_aii(pam_basal, DOSE_MAX, Emax, EC50, r)

# Titulação do ATHOS-3: início em 20, ajuste a cada 5 min em passos de até 15, máximo 80,
# alvo PAM ≥ 75, noradrenalina de fundo fixa.
# `ganho` (ng/kg/min por mmHg de erro) define o tamanho do passo; o documento só fixa o teto de 15.
# PK de um compartimento em dose equivalente x (= CL × concentração): dx/dt = k (u − x), k = ln2 / t½
function titular(p::Planta, Emax, r; ganho, dt = 1.0, alvo = PAM_ALVO_ATHOS)
    EC50  = ec50_aii(Emax)
    basal = pam_regime(p)
    k     = log(2) / (T_MEIA_AII * 60)
    decai = exp(-k * dt)
    n     = round(Int, HORIZONTE * 60 / dt)
    najuste = round(Int, INTERVALO_AJUSTE * 60 / dt)
    t = Vector{Float64}(undef, n + 1);  dose = similar(t);  pam = similar(t);  xs = similar(t)
    u = DOSE_INICIAL;  x = 0.0
    for i in 0:n
        if i > 0
            x = u + (x - u) * decai            # solução exata da PK com dose constante no passo
            if i % najuste == 0
                u = clamp(u + clamp(ganho * (alvo - pam_com_aii(basal, x, Emax, EC50, r)), -PASSO_MAX, PASSO_MAX),
                          0.0, DOSE_MAX)
            end
        end
        t[i + 1] = i * dt / 60;  dose[i + 1] = u;  xs[i + 1] = x;  pam[i + 1] = pam_com_aii(basal, x, Emax, EC50, r)
    end
    em(tmin) = pam[round(Int, tmin * 60 / dt) + 1]
    i75  = findfirst(>=(alvo - 0.5), pam)       # PAM que arredonda para o alvo
    isat = findfirst(>=(DOSE_MAX), dose)
    return (; t, dose, pam, x = xs, basal, EC50, r,
            dose_media = sum(dose[1:end-1]) / n, dose_final = dose[end],
            pam_1h = em(60), pam_3h = pam[end],
            t_alvo = i75 === nothing ? nothing : t[i75],
            t_saturacao = isat === nothing ? nothing : t[isat],
            F15 = (em(15) - basal) / (alvo - basal),
            pam_max = pam_max(basal, Emax, EC50, r),
            rp_final = p.Rp * pam[end] / basal)
end

# ------------------------------------------------------------------ noradrenalina + Ang II (troca de efeito)
#
# Do documento: dose de noradrenalina de 0 a 0,5 µg/kg/min (escolha de projeto), alta dose acima de 0,2
# (ATHOS-3), t½ de 2 a 6,8 min piorando com a gravidade (Beloeil 2005), PD de Hill sobre R_p somada à da
# Ang II, e o critério 4: no choque refratário a noradrenalina acima de 0,2 não atinge a PAM alvo.
# Premissas (não estão no documento, ficam em sliders no aplicativo):
#   * Emax da noradrenalina = fração do teto permitido pelo critério 4;
#   * EC50 da noradrenalina;
#   * regra de titulação da nora: abaixo de 65 sobe no passo máximo (15 × 0,5/80 por ajuste);
#     acima de 65 desce pela mesma regra proporcional da Ang II, escalada pela dose máxima; a nora persegue PAM 65 (alvo da SSC) e a Ang II persegue 75 (ATHOS-3), então entre
#     65 e 75 a Ang II sobe enquanto a nora desce;
#   * fase inicial de 60 min só com noradrenalina, partindo de zero, antes de iniciar a Ang II.

const DOSE_NE_MAX  = 0.5      # µg/kg/min
const DOSE_NE_ALTA = 0.2      # µg/kg/min (limiar de alta dose do ATHOS-3)
const T_MEIA_NE_FAIXA = (2.0, 6.8)   # min
const FASE_NE      = 60.0     # min só com noradrenalina antes da Ang II

# Efeito sobre R_p que a noradrenalina precisaria para levar a PAM a 65; o critério 4 exige Emax abaixo disso
teto_ne(p::Planta) = PAM_ALVO_SSC / pam_regime(p) - 1

function troca(p::Planta, Emax, r; ganho, frac_emax_ne, ec50_ne, t_meia_ne, dt = 1.0)
    EC50  = ec50_aii(Emax)
    EmaxN = frac_emax_ne * teto_ne(p)
    basal = pam_regime(p)
    escala = DOSE_NE_MAX / DOSE_MAX                       # converte a regra da Ang II para a nora
    decA  = exp(-log(2) / (T_MEIA_AII * 60) * dt)
    decN  = exp(-log(2) / (t_meia_ne * 60) * dt)
    n     = round(Int, (FASE_NE + HORIZONTE) * 60 / dt)
    najuste = round(Int, INTERVALO_AJUSTE * 60 / dt)
    t = Vector{Float64}(undef, n + 1)
    uN = similar(t); uA = similar(t); eN = similar(t); eA = similar(t); pam = similar(t)
    xNs = similar(t); xAs = similar(t)       # doses efetivas (PK), não as doses infundidas
    dN = 0.0; xN = 0.0; dA = 0.0; xA = 0.0
    pamf(xN, xA) = basal * (1 + hill(xN, EmaxN, ec50_ne) + hill(xA, Emax, r * EC50))
    for i in 0:n
        tm = i * dt / 60                     # min desde o início da noradrenalina
        if i > 0
            xN = dN + (xN - dN) * decN
            xA = dA + (xA - dA) * decA
            if i % najuste == 0
                P  = pamf(xN, xA)
                passoN = P < PAM_ALVO_SSC ? PASSO_MAX * escala :
                         -min(ganho * escala * (P - PAM_ALVO_SSC), PASSO_MAX * escala)
                dN = clamp(dN + passoN, 0.0, DOSE_NE_MAX)
                if tm >= FASE_NE
                    dA = clamp(dA + clamp(ganho * (PAM_ALVO_ATHOS - P), -PASSO_MAX, PASSO_MAX), 0.0, DOSE_MAX)
                end
            end
            if tm >= FASE_NE && dA == 0.0 && (i - 1) * dt / 60 < FASE_NE
                dA = DOSE_INICIAL                          # início da Ang II em 20 ng/kg/min
            end
        end
        t[i+1] = tm;  uN[i+1] = dN;  uA[i+1] = dA;  xNs[i+1] = xN;  xAs[i+1] = xA
        eN[i+1] = hill(xN, EmaxN, ec50_ne);  eA[i+1] = hill(xA, Emax, r * EC50);  pam[i+1] = pamf(xN, xA)
    end
    i0 = round(Int, FASE_NE * 60 / dt) + 1
    return (; t, nora = uN, aii = uA, x_nora = xNs, x_aii = xAs, efeito_nora = eN, efeito_aii = eA, pam, basal, EmaxN,
            pam_plato = pam[i0], nora_final = uN[end], aii_final = uA[end], pam_final = pam[end])
end
