# Aba 3 — Modelo Cardiovascular 3D & Circuito Equivalente Windkessel (WK3)
# Visualização esquemática 3D guiada pelas variáveis reais de modelo.jl + inspetor físico-matemático.

const GB3D = GLMakie.GeometryBasics
const Face3D = GLMakie.GeometryBasics.GLTriangleFace

function malha_tubo(pts::Vector{Point3f}, raios::Vector{Float32}; nseg::Int = 24)
    m = length(pts)
    verts = Vector{Point3f}(undef, m * nseg)
    faces = Face3D[]
    sizehint!(faces, 2 * (m - 1) * nseg)

    for i in 1:m
        tang = if i == 1
            pts[2] - pts[1]
        elseif i == m
            pts[m] - pts[m - 1]
        else
            pts[i + 1] - pts[i - 1]
        end
        nt = sqrt(tang[1]^2 + tang[2]^2 + tang[3]^2)
        t = nt > 1f-6 ? tang / nt : Point3f(0, 0, 1)

        ref = abs(t[2]) < 0.9f0 ? Point3f(0, 1, 0) : Point3f(1, 0, 0)
        n1 = Point3f(t[2] * ref[3] - t[3] * ref[2],
                     t[3] * ref[1] - t[1] * ref[3],
                     t[1] * ref[2] - t[2] * ref[1])
        nn1 = sqrt(n1[1]^2 + n1[2]^2 + n1[3]^2)
        n1 = nn1 > 1f-6 ? n1 / nn1 : Point3f(1, 0, 0)
        n2 = Point3f(t[2] * n1[3] - t[3] * n1[2],
                     t[3] * n1[1] - t[1] * n1[3],
                     t[1] * n1[2] - t[2] * n1[1])

        r = raios[i]
        for j in 1:nseg
            ang = Float32(2pi * (j - 1) / nseg)
            verts[(i - 1) * nseg + j] = pts[i] + r * (cos(ang) * n1 + sin(ang) * n2)
        end
    end

    for i in 1:(m - 1), j in 1:nseg
        jn = j == nseg ? 1 : j + 1
        a = (i - 1) * nseg + j
        b = (i - 1) * nseg + jn
        c = i * nseg + j
        d = i * nseg + jn
        push!(faces, Face3D(a, b, d))
        push!(faces, Face3D(a, d, c))
    end
    return GB3D.Mesh(verts, faces)
end

function malha_elipsoide(centro::Point3f, rx::Float32, ry::Float32, rz::Float32; nu::Int = 28, nv::Int = 20)
    verts = Point3f[]
    faces = Face3D[]
    for iv in 0:nv
        phi = Float32(pi * iv / nv)
        for iu in 0:(nu - 1)
            theta = Float32(2pi * iu / nu)
            x = centro[1] + rx * sin(phi) * cos(theta)
            y = centro[2] + ry * sin(phi) * sin(theta)
            z = centro[3] + rz * cos(phi)
            push!(verts, Point3f(x, y, z))
        end
    end
    for iv in 0:(nv - 1), iu in 0:(nu - 1)
        iun = iu == nu - 1 ? 0 : iu + 1
        a = iv * nu + iu + 1
        b = iv * nu + iun + 1
        c = (iv + 1) * nu + iu + 1
        d = (iv + 1) * nu + iun + 1
        push!(faces, Face3D(a, b, d))
        push!(faces, Face3D(a, d, c))
    end
    return GB3D.Mesh(verts, faces)
end

function info_componente_wk(comp::Int, d, inst)
    k = inst.k
    p = d.plantas[k]
    s = serie(d, k)
    i = indice(s, inst.tm)
    e_tot = s.e_aii[i] + (s.e_ne === nothing ? 0.0 : s.e_ne[i])
    rp_atual = p.Rp * (1.0 + e_tot)
    zc_atual = zc(rp_atual)
    dc_atual = dc_mls(p)
    b = inst.bat[k]
    pc_med = rp_atual * dc_atual + P_VEN

    if comp == 1
        return (
            nome = "1 · Coração (Bomba Ventricular Esquerda)",
            resumo = "Qᵢₙ(t)  ·  FC = " * f0(p.FC) * " bpm  ·  VS = " * f1(p.VS) * " ml  ·  DC = " * f2(dc_lmin(p)) * " L/min",
            eq = "Qᵢₙ(t) = (π · VS / 2 · Tₑⱼ) · sin(π · t / Tₑⱼ)    para 0 ≤ t < Tₑⱼ",
            papel = "Gera o fluxo pulsátil de ejeção durante a sístole (Tₑⱼ = 0,413 − 0,0017·FC) e permanece nulo na diástole. No choque séptico, a taquicardia compensatória eleva o débito cardíaco frente à vasoplegia.",
            influencia = "Define o volume ejetado por batimento. Multiplicado pela resistência total (Rₚ + Z꜀), determina a pressão média: PAM = (Rₚ + Z꜀) · DC + Pᵥₑₙ."
        )
    elseif comp == 2
        return (
            nome = "2 · Impedância Característica da Aorta (Z꜀)",
            resumo = "Z꜀ = 0,06 · Rₚ  ·  Valor atual: " * replace(@sprintf("%.4f", zc_atual), "." => ",") * " mmHg·s/ml",
            eq = "Z꜀ = 0,06 · Rₚ      |      Pₐₒ(t) = P꜀(t) + Z꜀ · Qᵢₙ(t)",
            papel = "Representa a oposição inercial e elástica imediata da raiz da aorta ao fluxo rápido ejetado pelo ventrículo esquerdo no início da sístole.",
            influencia = "Adiciona a queda resistiva instantânea Z꜀ · Qᵢₙ(t) sobre P꜀(t) durante a sístole, elevando o pico sistólico (PAS) e contribuindo levemente para a PAM."
        )
    elseif comp == 3
        return (
            nome = "3 · Complacência Arterial (C)",
            resumo = "C = " * f2(p.C) * " ml/mmHg  ·  Fator idade: " * f2(PARAM.fcx[]) * "×",
            eq = "C = ΔV / ΔP      |      C · (dP꜀/dt) = Qᵢₙ(t) − (P꜀ − Pᵥₑₙ) / Rₚ",
            papel = "Representa a capacidade elástica das grandes artérias de armazenar volume na sístole e devolvê-lo na diástole, amortecendo a pulsatilidade do fluxo.",
            influencia = "Não altera a PAM média em regime, mas controla a pressão de pulso (PP = PAS − PAD). No idoso (C menor), a PAS sobe e a PAD cai (PP atual: " * f0(b.PP) * " mmHg)."
        )
    elseif comp == 4
        return (
            nome = "4 · Resistência Vascular Periférica (Rₚ)",
            resumo = "Rₚ atual = " * replace(@sprintf("%.3f", rp_atual), "." => ",") * " mmHg·s/ml  ·  Basal: " * replace(@sprintf("%.3f", p.Rp), "." => ","),
            eq = "Rₚ(t) = Rₚ,₀ · (1 + E_NE + E_AII)      |      E(x) = Eₘₐₓ · x / (r · EC₅₀ + x)",
            papel = "Representa o tônus das arteríolas sistêmicas, principal alvo farmacológico da Noradrenalina (receptores α₁) e da Angiotensina II (receptores AT₁).",
            influencia = "Principal determinante da PAM. No choque vasodilatador, Rₚ colapsa; o aumento de Rₚ pelos vasopressores eleva proporcionalmente PAM, PAS e PAD."
        )
    elseif comp == 5
        return (
            nome = "5 · Pressão Aórtica Instantânea (Pₐₒ)",
            resumo = "PAM = " * f1(inst.pam) * " mmHg  ·  PAS/PAD = " * f0(b.PAS) * "/" * f0(b.PAD) * " mmHg  ·  PP = " * f0(b.PP) * " mmHg",
            eq = "Pₐₒ(t) = P꜀(t) + Z꜀ · Qᵢₙ(t)",
            papel = "Pressão arterial na raiz da aorta ao longo do ciclo cardíaco, resultante da soma entre a pressão do reservatório elástico (P꜀) e a componente proximal (Z꜀ · Qᵢₙ).",
            influencia = "Variável clínica monitorada. Meta terapêutica do estudo ATHOS-3: PAM ≥ 75 mmHg (referência Surviving Sepsis Campaign: PAM ≥ 65 mmHg)."
        )
    else
        return (
            nome = "6 · Reservatório Elástico & Referência Venosa (P꜀ / Pᵥₑₙ)",
            resumo = "P꜀ médio ≈ " * f1(pc_med) * " mmHg  ·  Pᵥₑₙ = " * f0(P_VEN) * " mmHg  ·  τ = " * f2(rp_atual * p.C) * " s",
            eq = "dP꜀/dt = [Qᵢₙ(t) − (P꜀ − Pᵥₑₙ) / Rₚ] / C    (Euler dt = 0,1 ms)",
            papel = "P꜀(t) é a variável de estado da EDO do Windkessel (pressão no reservatório arterial). Pᵥₑₙ = 0 mmHg é a pressão venosa central de referência adotada.",
            influencia = "Na diástole (Qᵢₙ = 0), P꜀ decai exponencialmente com constante de tempo τ = Rₚ · C, sustentando a perfusão diastólica e coronariana."
        )
    end
end

function aba_modelo3d(corpo)
    SEL_COMP = Observable(3)
    NL = string('\n')

    # ================================================================
    # PAINEL ESQUERDO (350 px, limpo e bem espaçado)
    # ================================================================
    painel = corpo[1, 1] = GridLayout(tellheight = false, valign = :top)
    Label(painel[1, 1], "Modelo Cardiovascular 3D"; fontsize = 18, font = :bold, halign = :left)
    Label(painel[2, 1], lift((p, u) -> (p == 1 ? "Adulto (46 anos)" : "Idoso (70 anos)") *
                                       "  ·  " * ["Sem IECA/BRA", "IECA", "BRA"][u], PAC, USO);
          halign = :left, color = AZUL, font = :bold, fontsize = 13)

    Label(painel[3, 1], "Elemento inspecionado"; halign = :left, font = :bold, fontsize = 13, padding = (0, 0, 0, 6))
    comp_grid = painel[4, 1] = GridLayout()
    nomes_comp = [
        "1 · Coração (Qᵢₙ)",
        "2 · Impedância (Z꜀)",
        "3 · Complacência (C)",
        "4 · Resistência (Rₚ)",
        "5 · Pressão (Pₐₒ)",
        "6 · Reservatório (P꜀)"
    ]
    bts_comp = Any[]
    for idx in 1:6
        r = div(idx - 1, 2) + 1
        c = mod(idx - 1, 2) + 1
        bt = Button(comp_grid[r, c];
                    label = nomes_comp[idx],
                    width = 165,
                    fontsize = 12,
                    buttoncolor = idx == SEL_COMP[] ? ATIVA : INATIVA,
                    labelcolor = idx == SEL_COMP[] ? :black : :white)
        on(_ -> (SEL_COMP[] = idx), bt.clicks)
        push!(bts_comp, bt)
    end
    on(SEL_COMP) do s
        for (idx, bt) in enumerate(bts_comp)
            bt.buttoncolor[] = idx == s ? ATIVA : INATIVA
            bt.labelcolor[]  = idx == s ? :black : :white
        end
    end
    rowgap!(comp_grid, 6)
    colgap!(comp_grid, 8)

    Label(painel[5, 1], "Parâmetros vasculares"; halign = :left, font = :bold, fontsize = 13, padding = (0, 0, 0, 8))
    sg_fis = SliderGrid(painel[6, 1],
        (label = "Rₚ com a idade (×)", range = 1.10:0.01:1.20, startvalue = PARAM.frp[], format = "{:.2f}"),
        (label = "C com a idade (×)", range = 0.25:0.01:0.50, startvalue = PARAM.fcx[], format = "{:.2f}"),
        (label = "FC no choque (bpm)", range = 63:1:120, startvalue = PARAM.fc[], format = "{:.0f}"))
    for (s, o) in zip(sg_fis.sliders, (PARAM.frp, PARAM.fcx, PARAM.fc))
        on(s.value) do v
            o[] != v && (o[] = v)
        end
        on(o) do v
            s.value[] != v && set_close_to!(s, v)
        end
    end

    Label(painel[7, 1], "Leitura visual"; halign = :left, font = :bold, fontsize = 13, padding = (0, 0, 0, 8))
    Label(painel[8, 1], join([
        "• Coração pulsa conforme Qᵢₙ(t) e VS",
        "• Câmara central expande com C e P꜀(t)",
        "• Arteríolas fecham conforme Rₚ(t) sobe",
        "• Arraste sobre o 3D para girar ou dar zoom"
    ], NL); halign = :left, justification = :left, color = GRIS, fontsize = 12, tellwidth = false, word_wrap = true)

    Label(painel[9, 1], "Representação esquemática 3D guiada em tempo real pelas equações 0D do modelo (não é CFD).";
          halign = :left, justification = :left, color = GRIS, fontsize = 11, tellwidth = false, word_wrap = true,
          padding = (0, 0, 0, 6))
    rowgap!(painel, 6)

    # ================================================================
    # ÁREA PRINCIPAL (DIREITA)
    # ================================================================
    direita = corpo[1, 2] = GridLayout()
    topo = direita[1, 1] = GridLayout()

    # ----------------------------------------------------------------
    # A. CENA ANATÔMICA 3D
    # ----------------------------------------------------------------
    titulo_3d = Observable("A · Sistema Cardiovascular Esquemático 3D")
    ax3 = Axis3(topo[1, 1];
                title = titulo_3d,
                titlesize = 14,
                aspect = :data,
                azimuth = 1.35pi,
                elevation = 0.14pi,
                perspectiveness = 0.35)
    hidedecorations!(ax3)
    hidespines!(ax3)

    mesh_coracao = Observable(malha_elipsoide(Point3f(-2.6, 0.0, 0.0), 0.75f0, 0.65f0, 0.85f0))
    cor_coracao  = Observable(VERMELHO)

    pts_aorta = [
        Point3f(-2.05, 0.0, 0.45),
        Point3f(-1.85, 0.0, 1.15),
        Point3f(-1.35, 0.0, 1.65),
        Point3f(-0.65, 0.0, 1.75),
        Point3f( 0.00, 0.0, 1.45)
    ]
    mesh_aorta = Observable(malha_tubo(pts_aorta, fill(0.24f0, length(pts_aorta))))
    cor_aorta  = Observable(LARANJA)

    pts_comp = [
        Point3f(0.00, 0.0, 1.45),
        Point3f(0.25, 0.0, 0.90),
        Point3f(0.35, 0.0, 0.20),
        Point3f(0.35, 0.0, -0.50),
        Point3f(0.25, 0.0, -1.05)
    ]
    mesh_comp = Observable(malha_tubo(pts_comp, fill(0.34f0, length(pts_comp))))
    cor_comp  = Observable(AZUL)

    pts_ramo_sup = [Point3f(-1.10, 0.0, 1.70), Point3f(-1.05, 0.25, 2.25), Point3f(-0.95, 0.40, 2.65)]
    pts_ramo_esq = [Point3f(0.25, 0.0, -1.05), Point3f(0.75, -0.45, -1.55), Point3f(1.35, -0.75, -1.85)]
    pts_ramo_dir = [Point3f(0.25, 0.0, -1.05), Point3f(0.75,  0.45, -1.55), Point3f(1.35,  0.75, -1.85)]
    mesh_ga1 = malha_tubo(pts_ramo_sup, Float32[0.14, 0.12, 0.10])
    mesh_ga2 = malha_tubo(pts_ramo_esq, Float32[0.18, 0.15, 0.13])
    mesh_ga3 = malha_tubo(pts_ramo_dir, Float32[0.18, 0.15, 0.13])

    rotas_perif = [
        [Point3f(-0.95, 0.40, 2.65), Point3f(-0.40, 0.65, 3.05), Point3f(0.35, 0.80, 3.20)],
        [Point3f(1.35, -0.75, -1.85), Point3f(2.05, -1.05, -2.05), Point3f(2.75, -1.20, -2.10)],
        [Point3f(1.35, -0.75, -1.85), Point3f(2.10, -0.45, -1.55), Point3f(2.80, -0.35, -1.40)],
        [Point3f(1.35,  0.75, -1.85), Point3f(2.10,  0.45, -1.55), Point3f(2.80,  0.35, -1.40)],
        [Point3f(1.35,  0.75, -1.85), Point3f(2.05,  1.05, -2.05), Point3f(2.75,  1.20, -2.10)]
    ]
    meshes_perif = [Observable(malha_tubo(r, Float32[0.10, 0.08, 0.06])) for r in rotas_perif]
    cor_perif    = Observable(VERDE)

    mesh!(ax3, mesh_coracao; color = cor_coracao, shading = FastShading)
    mesh!(ax3, mesh_aorta; color = cor_aorta, shading = FastShading)
    mesh!(ax3, mesh_comp; color = cor_comp, shading = FastShading, transparency = true)
    mesh!(ax3, mesh_ga1; color = RGBf(0.75, 0.35, 0.38), shading = FastShading)
    mesh!(ax3, mesh_ga2; color = RGBf(0.75, 0.35, 0.38), shading = FastShading)
    mesh!(ax3, mesh_ga3; color = RGBf(0.75, 0.35, 0.38), shading = FastShading)
    for mp in meshes_perif
        mesh!(ax3, mp; color = cor_perif, shading = FastShading)
    end

    text!(ax3, Point3f(-2.6, 0.0, -1.15); text = "Coração (Qᵢₙ)", color = TEXTO, fontsize = 12, align = (:center, :top))
    text!(ax3, Point3f(-1.25, 0.0, 2.05); text = "Aorta (Z꜀)", color = LARANJA, fontsize = 12, align = (:center, :bottom))
    text!(ax3, Point3f(0.95, 0.0, 0.25); text = "Complacência (C)", color = AZUL, fontsize = 12, align = (:left, :center))
    text!(ax3, Point3f(2.85, 0.0, -1.75); text = "Leito Periférico (Rₚ)", color = VERDE, fontsize = 12, align = (:left, :center))
    limits!(ax3, -3.7, 4.1, -2.0, 2.0, -2.5, 3.6)

    # ----------------------------------------------------------------
    # B. CIRCUITO EQUIVALENTE WINDKESSEL 3 ELEMENTOS (WK3) — SEM SOBREPOSIÇÃO
    # ----------------------------------------------------------------
    circ_card = topo[1, 2] = GridLayout()
    try
        Box(circ_card[1, 1]; color = PAINEL, strokecolor = BORDA, strokewidth = 1.2, cornerradius = 10)
    catch
        Box(circ_card[1, 1]; color = PAINEL, strokecolor = BORDA, strokewidth = 1.2)
    end

    ax_circ = Axis(circ_card[1, 1];
                   title = "B · Circuito Equivalente Windkessel (WK3)",
                   titlesize = 14,
                   subtitle = "Valores instantâneos do paciente e cenário selecionados",
                   subtitlesize = 12,
                   backgroundcolor = :transparent,
                   limits = (-2.0, 11.6, 0.2, 6.8))
    hidedecorations!(ax_circ)
    hidespines!(ax_circ)

    # Fios do circuito (malha esquerda em x = 2.2, capacitor em x = 6.2, resistor R_p em x = 9.0)
    fio = RGBf(0.55, 0.63, 0.74)
    lines!(ax_circ, [2.2, 9.0], [1.3, 1.3]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [2.2, 2.2], [1.3, 2.55]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [2.2, 2.2], [3.85, 5.1]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [2.2, 3.3], [5.1, 5.1]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [5.1, 6.2], [5.1, 5.1]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [6.2, 6.2], [5.1, 3.65]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [6.2, 6.2], [2.95, 1.3]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [6.2, 9.0], [5.1, 5.1]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [9.0, 9.0], [5.1, 4.25]; color = fio, linewidth = 2.5)
    lines!(ax_circ, [9.0, 9.0], [2.15, 1.3]; color = fio, linewidth = 2.5)

    # Símbolo 1: Fonte de fluxo Q_in
    ang_c = range(0, 2pi, length = 60)
    lines!(ax_circ, 2.2 .+ 0.62 .* cos.(ang_c), 3.2 .+ 0.62 .* sin.(ang_c); color = VERMELHO, linewidth = 3)
    arrows!(ax_circ, [Point2f(2.2, 2.8)], [Vec2f(0.0, 0.75)]; color = VERMELHO, linewidth = 2.5, arrowsize = 13)

    # Símbolo 2: Resistor Z_c
    zx = [3.3, 3.5, 3.75, 4.05, 4.35, 4.65, 4.9, 5.1]
    zy = [5.1, 5.4, 4.8, 5.4, 4.8, 5.4, 4.8, 5.1]
    lines!(ax_circ, zx, zy; color = LARANJA, linewidth = 3)

    # Símbolo 3: Capacitor C
    lines!(ax_circ, [5.6, 6.8], [3.65, 3.65]; color = AZUL, linewidth = 3.5)
    lines!(ax_circ, [5.6, 6.8], [2.95, 2.95]; color = AZUL, linewidth = 3.5)

    # Símbolo 4: Resistor R_p
    rx = [9.0, 9.3, 8.7, 9.3, 8.7, 9.3, 8.7, 9.0]
    ry = [4.25, 4.05, 3.75, 3.4, 3.05, 2.7, 2.35, 2.15]
    lines!(ax_circ, rx, ry; color = VERDE, linewidth = 3)

    # Nós P_ao e P_c
    scatter!(ax_circ, [2.75, 6.2], [5.1, 5.1]; color = [LARANJA, AZUL], markersize = 12)

    txt_qin = lift((i, d) -> join([
        "Qᵢₙ (Coração)",
        "DC = " * f2(dc_lmin(d.plantas[i.k])) * " L/min",
        "VS = " * f1(d.plantas[i.k].VS) * " ml"
    ], NL), INST, DADOS)
    txt_pao = lift(i -> "Nó Pₐₒ = " * f1(i.pam) * " mmHg" * NL * "(" * f0(i.bat[i.k].PAS) * "/" * f0(i.bat[i.k].PAD) * ")", INST)
    txt_zc  = lift((i, d) -> begin
        s = serie(d, i.k); idx = indice(s, i.tm)
        et = s.e_aii[idx] + (s.e_ne === nothing ? 0.0 : s.e_ne[idx])
        rp = d.plantas[i.k].Rp * (1.0 + et)
        "Z꜀ = " * replace(@sprintf("%.4f", zc(rp)), "." => ",") * NL * "mmHg·s/ml"
    end, INST, DADOS)
    txt_pc  = lift((i, d) -> begin
        s = serie(d, i.k); idx = indice(s, i.tm)
        et = s.e_aii[idx] + (s.e_ne === nothing ? 0.0 : s.e_ne[idx])
        rp = d.plantas[i.k].Rp * (1.0 + et)
        pcm = rp * dc_mls(d.plantas[i.k]) + P_VEN
        "Nó P꜀ ≈ " * f1(pcm) * " mmHg"
    end, INST, DADOS)
    txt_c   = lift((i, d) -> "C = " * f2(d.plantas[i.k].C) * NL * "ml/mmHg", INST, DADOS)
    txt_rp  = lift((i, d) -> begin
        s = serie(d, i.k); idx = indice(s, i.tm)
        et = s.e_aii[idx] + (s.e_ne === nothing ? 0.0 : s.e_ne[idx])
        rp = d.plantas[i.k].Rp * (1.0 + et)
        "Rₚ = " * replace(@sprintf("%.3f", rp), "." => ",") * NL * "mmHg·s/ml"
    end, INST, DADOS)

    # Posições livres de qualquer linha ou símbolo
    text!(ax_circ, 1.35, 3.20; text = txt_qin, color = VERMELHO, fontsize = 11, align = (:right, :center))
    text!(ax_circ, 2.75, 5.55; text = txt_pao, color = LARANJA, fontsize = 11, font = :bold, align = (:center, :bottom))
    text!(ax_circ, 4.20, 4.55; text = txt_zc, color = LARANJA, fontsize = 11, align = (:center, :top))
    text!(ax_circ, 6.20, 5.55; text = txt_pc, color = AZUL, fontsize = 11, font = :bold, align = (:center, :bottom))
    text!(ax_circ, 6.95, 3.30; text = txt_c, color = AZUL, fontsize = 11, align = (:left, :center))
    text!(ax_circ, 9.45, 3.20; text = txt_rp, color = VERDE, fontsize = 11, align = (:left, :center))
    text!(ax_circ, 5.60, 0.90; text = "Referência Venosa: Pᵥₑₙ = 0,0 mmHg (Melo & Peters, 1999)",
          color = GRIS, fontsize = 11, align = (:center, :top))

    colsize!(topo, 1, Relative(0.48))
    colsize!(topo, 2, Relative(0.52))
    colgap!(topo, 16)

    # ================================================================
    # LINHA INFERIOR: [1, 1] CARTÃO INSPETOR  |  [1, 2] ONDA P_ao(t)
    # ================================================================
    base = direita[2, 1] = GridLayout()

    insp_box = base[1, 1] = GridLayout()
    try
        Box(insp_box[1:5, 1]; color = CARTAO, strokecolor = BORDA, strokewidth = 1.5, cornerradius = 10)
    catch
        Box(insp_box[1:5, 1]; color = CARTAO, strokecolor = BORDA, strokewidth = 1.5)
    end
    info_obs = lift((c, d, inst) -> info_componente_wk(c, d, inst), SEL_COMP, DADOS, INST)

    # tellwidth = false + word_wrap = true impede o texto de invadir o gráfico da direita
    Label(insp_box[1, 1], lift(x -> x.nome, info_obs);
          fontsize = 15, font = :bold, color = AZUL, halign = :left, padding = (16, 16, 10, 2),
          tellwidth = false, word_wrap = true)
    Label(insp_box[2, 1], lift(x -> x.resumo, info_obs);
          fontsize = 12, font = :bold, color = TEXTO, halign = :left, padding = (16, 16, 2, 2),
          tellwidth = false, word_wrap = true)
    Label(insp_box[3, 1], lift(x -> "Equação:  " * x.eq, info_obs);
          fontsize = 12, font = :bold, color = AMARELO, halign = :left, padding = (16, 16, 2, 4),
          tellwidth = false, word_wrap = true)
    Label(insp_box[4, 1], lift(x -> "Papel fisiológico: " * x.papel, info_obs);
          fontsize = 12, color = TEXTO, halign = :left, justification = :left, word_wrap = true,
          padding = (16, 16, 2, 2), tellwidth = false)
    Label(insp_box[5, 1], lift(x -> "Impacto hemodinâmico: " * x.influencia, info_obs);
          fontsize = 12, color = VERDE, halign = :left, justification = :left, word_wrap = true,
          padding = (16, 16, 2, 10), tellwidth = false)
    rowgap!(insp_box, 3)

    sub_onda3d = Observable("")
    ax_onda3d = Axis(base[1, 2];
                     title = "C · Onda Pulsátil Pₐₒ(t)",
                     titlesize = 14,
                     subtitle = sub_onda3d,
                     subtitlesize = 12,
                     xlabel = "tempo no ciclo (s)",
                     xlabelsize = 12,
                     ylabel = "Pₐₒ (mmHg)",
                     ylabelsize = 12)
    obs_onda3d  = Observable(Point2f[(0, 0)])
    obs_marca3d = Observable([Point2f(0, 0)])
    lines!(ax_onda3d, obs_onda3d; color = LARANJA, linewidth = 2.8)
    scatter!(ax_onda3d, obs_marca3d; color = :white, strokecolor = VERMELHO, strokewidth = 2.5, markersize = 12)

    colsize!(base, 1, Relative(0.60))
    colsize!(base, 2, Relative(0.40))
    colgap!(base, 18)

    rowsize!(direita, 1, Relative(0.60))
    rowsize!(direita, 2, Relative(0.40))
    rowgap!(direita, 16)

    # ================================================================
    # ATUALIZAÇÃO DINÂMICA DAS GEOMETRIAS 3D E DA ONDA
    # ================================================================
    function atualizar_geometrias_3d!()
        inst = INST[]; d = DADOS[]
        (inst === nothing || d === nothing) && return
        k = inst.k
        p = d.plantas[k]
        s = serie(d, k)
        idx = indice(s, inst.tm)
        et = s.e_aii[idx] + (s.e_ne === nothing ? 0.0 : s.e_ne[idx])
        rp_atual = p.Rp * (1.0 + et)

        T = 60.0 / p.FC
        Tej = t_ejecao(p.FC)
        tc = mod(TBEAT[], T)
        q_norm = clamp(q_in(tc, p.VS, Tej) / max(1.0, pi * p.VS / (2 * Tej)), 0.0, 1.0)

        fator_vs = Float32(clamp(p.VS / 100.0, 0.75, 1.35))
        pulso_c = Float32(1.0 - 0.16 * q_norm) * sqrt(fator_vs)
        mesh_coracao[] = malha_elipsoide(Point3f(-2.6, 0.0, 0.0), 0.78f0 * pulso_c, 0.68f0 * pulso_c, 0.90f0 * pulso_c)

        fator_c = Float32(clamp(sqrt(p.C / C_ADULTO), 0.45, 1.15))
        exp_c   = Float32(1.0 + 0.14 * q_norm * (p.C / C_ADULTO))
        r_comp  = 0.42f0 * fator_c * exp_c
        mesh_comp[] = malha_tubo(pts_comp, fill(r_comp, length(pts_comp)))

        fator_rp = Float32(clamp((RP_CHOQUE_ADULTO / rp_atual)^0.35, 0.45, 1.25))
        for (j, rota) in enumerate(rotas_perif)
            meshes_perif[j][] = malha_tubo(rota, Float32[0.11, 0.085, 0.065] .* fator_rp)
        end

        sc = SEL_COMP[]
        cor_coracao[] = sc == 1 ? AMARELO : VERMELHO
        cor_aorta[]   = sc in (2, 5) ? AMARELO : LARANJA
        cor_comp[]    = sc in (3, 6) ? AMARELO : comalfa(AZUL, 0.85)
        cor_perif[]   = sc == 4 ? AMARELO : VERDE
        return
    end

    function atualizar_onda_3d!()
        inst = INST[]; d = DADOS[]
        (inst === nothing || d === nothing) && return
        k = inst.k
        b = inst.bat[k]
        T = 60.0 / d.plantas[k].FC
        obs_onda3d[] = vcat(Point2f.(b.t, b.pao), Point2f.(b.t .+ T, b.pao))
        limits!(ax_onda3d, 0, 2 * T, 0, max(140.0, b.PAS * 1.25))
        sub_onda3d[] = ROTULOS[k] * " · " * f0(b.PAS) * "/" * f0(b.PAD) * " mmHg (PP " * f0(b.PP) * ")"
        titulo_3d[]  = "A · Sistema Cardiovascular Esquemático 3D" * NL *
                       ROTULOS[k] * " · t = " * f0(inst.tm) * " min · PAM " * f1(inst.pam) * " mmHg"
        return
    end

    function atualizar_marca_3d!()
        inst = INST[]; d = DADOS[]
        (inst === nothing || d === nothing) && return
        b = inst.bat[inst.k]
        T = 60.0 / d.plantas[inst.k].FC
        x = mod(TBEAT[], 2 * T)
        i = clamp(round(Int, mod(x, T) / DT_EULER) + 1, 1, length(b.pao))
        obs_marca3d[] = [Point2f(x, b.pao[i])]
        atualizar_geometrias_3d!()
        return
    end

    on(_ -> (atualizar_onda_3d!(); atualizar_geometrias_3d!()), INST)
    on(_ -> (atualizar_onda_3d!(); atualizar_geometrias_3d!()), DADOS)
    on(_ -> atualizar_geometrias_3d!(), SEL_COMP)
    on(_ -> atualizar_marca_3d!(), TBEAT)

    atualizar_onda_3d!()
    atualizar_geometrias_3d!()
    atualizar_marca_3d!()
    colsize!(corpo, 1, Fixed(350))
    return
end