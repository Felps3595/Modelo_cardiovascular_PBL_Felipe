# Aba 2 — Tratamento Farmacológico (PK/PD da Noradrenalina e Angiotensina II)

function aba_tratamento(corpo)
    painel = corpo[1, 1] = GridLayout(tellheight = false, valign = :top)
    Label(painel[1, 1], "Tratamento Farmacológico"; fontsize = 18, font = :bold, halign = :left)
    Label(painel[2, 1], lift((p, u) -> "Perfil ativo: " * (p == 1 ? "Adulto (46 a)" : "Idoso (70 a)") *
                                       " · " * ["Sem IECA/BRA", "IECA", "BRA"][u], PAC, USO);
          halign = :left, color = AZUL, font = :bold, fontsize = 13)

    lt = painel[3, 1] = GridLayout(halign = :left)
    tg = Toggle(lt[1, 1]; active = FONTE[] == :troca)
    Label(lt[1, 2], "simular desmame da noradrenalina\n(ativa série :troca nas séries temporais)";
          halign = :left, justification = :left, fontsize = 12)
    on(tg.active) do a
        FONTE[] != (a ? :troca : :titular) && (FONTE[] = a ? :troca : :titular)
    end
    on(FONTE) do f
        tg.active[] != (f == :troca) && (tg.active[] = (f == :troca))
    end

    Label(painel[4, 1], "Angiotensina II (ATHOS-3 / Leisman 2024)"; halign = :left, font = :bold, color = LARANJA)
    sg_aii = SliderGrid(painel[5, 1],
        (label = "Emax Ang II", range = 0.30:0.05:1.00, startvalue = PARAM.Emax[], format = "{:.2f}"),
        (label = "Ganho titulação", range = 0.3:0.1:2.0, startvalue = PARAM.ganho[], format = "{:.1f}"))
    on(sg_aii.sliders[1].value) do v
        PARAM.Emax[] != v && (PARAM.Emax[] = v)
    end
    on(PARAM.Emax) do v
        sg_aii.sliders[1].value[] != v && set_close_to!(sg_aii.sliders[1], v)
    end
    on(sg_aii.sliders[2].value) do v
        PARAM.ganho[] != v && (PARAM.ganho[] = v)
    end
    on(PARAM.ganho) do v
        sg_aii.sliders[2].value[] != v && set_close_to!(sg_aii.sliders[2], v)
    end

    Label(painel[6, 1], lift(d -> "EC50 calibrada: " * f1(d.EC50) * " ng/kg/min · t½: " * f1(T_MEIA_AII) * " min\n" *
                                  "Fator r (BRA): " * f1(d.rb) * "  →  EC50' = " * f1(d.rb * d.EC50) * " ng/kg/min", DADOS);
          halign = :left, justification = :left, color = GRIS, fontsize = 12)

    Label(painel[7, 1], "Noradrenalina (premissas do grupo)"; halign = :left, font = :bold, color = AZUL)
    sg_ne = SliderGrid(painel[8, 1],
        (label = "Emax NE (fração teto)", range = 0.20:0.05:0.95, startvalue = PARAM.ne_frac[], format = "{:.2f}"),
        (label = "EC50 NE (µg/kg/min)", range = 0.05:0.05:0.50, startvalue = PARAM.ne_ec50[], format = "{:.2f}"),
        (label = "t½ NE (min)", range = 2.0:0.2:6.8, startvalue = PARAM.ne_t12[], format = "{:.1f}"))
    for (s, o) in zip(sg_ne.sliders, (PARAM.ne_frac, PARAM.ne_ec50, PARAM.ne_t12))
        on(s.value) do v
            o[] != v && (o[] = v)
        end
    end

    Label(painel[9, 1], AVISO_NE; halign = :left, justification = :left, color = VERMELHO, fontsize = 11,
          word_wrap = true, width = 330)
    Label(painel[10, 1], "Nota: na simulação base (:titular), a noradrenalina é tratada como dose de fundo fixa " *
                        "já incorporada em R_p,choque. Ligue o desmame acima para observar a troca dinâmica (0–240 min).";
          halign = :left, justification = :left, color = GRIS, fontsize = 11, word_wrap = true, width = 330)
    rowgap!(painel, 6)

    direita = corpo[1, 2] = GridLayout()

    cards = direita[1, 1] = GridLayout()
    defs = [
        ("Dose NE",
         lift(i -> i.nora === nothing ? "fixa" : f2(i.nora), INST),
         lift(i -> i.nora === nothing ? ">0,2 µg/kg/min" : "µg/kg/min", INST),
         AZUL),
        ("Emax / EC50 NE",
         lift(d -> f2(d.trocas[K[]].EmaxN), DADOS),
         lift((ec, t12) -> "EC50 " * f2(ec) * " · t½ " * f1(t12) * "m", PARAM.ne_ec50, PARAM.ne_t12),
         AZUL),
        ("Efeito NE (ΔR_p)",
         lift(i -> i.e_ne === nothing ? "base" : "+" * f1(100 * i.e_ne) * "%", INST),
         lift((i, d) -> i.e_ne === nothing ? "R_p,0 = " * f2(d.plantas[i.k].Rp) :
                                             "+" * f2(d.plantas[i.k].Rp * i.e_ne) * " mmHg·s/ml", INST, DADOS),
         AZUL),
        ("Dose Ang II",
         lift(i -> f1(i.aii), INST),
         lift(d -> "ng/kg/min · t½ " * f1(T_MEIA_AII) * "m", DADOS),
         LARANJA),
        ("Fator r / EC50'",
         lift(d -> f1(d.rs[K[]]) * "×", DADOS),
         lift(d -> "EC50' " * f1(d.rs[K[]] * d.EC50) * " ng", DADOS),
         lift(k -> COR[k], K)),
        ("Efeito Ang II",
         lift(i -> "+" * f1(100 * i.e_aii) * "%", INST),
         lift((i, d) -> "ΔR_p +" * f2(d.plantas[i.k].Rp * i.e_aii) * " mmHg·s/ml", INST, DADOS),
         LARANJA)
    ]
    for (j, (tit, val, apoio, cor)) in enumerate(defs)
        sub = cards[1, j] = GridLayout()
        cartao!(sub, tit, val, apoio, cor)
    end
    colgap!(cards, 10)

    gg = direita[2, 1] = GridLayout()

    ax_hill_ne = Axis(gg[1, 1];
                      title = "A · Curva de Hill: Noradrenalina (premissa)",
                      titlesize = 13,
                      subtitle = "E_NE(x) = Emax_NE · x / (EC50_NE + x)",
                      subtitlesize = 11,
                      xlabel = "dose efetiva NE (µg/kg/min)",
                      xlabelsize = 12,
                      ylabel = "efeito relativo E_NE",
                      ylabelsize = 12)
    xs_ne = collect(Float32, range(0.0, DOSE_NE_MAX, length = 200))
    obs_hill_ne = Observable(Point2f[(0, 0)])
    obs_ec50_ne = Observable([Point2f(0, 0)])
    obs_pt_ne   = Observable([Point2f(NaN, NaN)])
    lines!(ax_hill_ne, obs_hill_ne; color = AZUL, linewidth = 3, label = "curva Hill NE")
    vlines!(ax_hill_ne, [DOSE_NE_ALTA]; color = GRIS, linestyle = :dash, linewidth = 1.2)
    scatter!(ax_hill_ne, obs_ec50_ne; color = AMARELO, markersize = 10, label = "EC50_NE")
    scatter!(ax_hill_ne, obs_pt_ne; color = :white, strokecolor = AZUL, strokewidth = 2.5, markersize = 12, label = "instante atual")
    axislegend(ax_hill_ne; position = :lt, framevisible = false, labelsize = 10, patchsize = (16, 10), rowgap = 1)

    ax_hill_aii = Axis(gg[1, 2];
                       title = "B · Curvas de Hill: Angiotensina II (deslocamento por BRA)",
                       titlesize = 13,
                       subtitle = "E_AII(x) = Emax · x / (r · EC50 + x)  |  EC50' = r × EC50",
                       subtitlesize = 11,
                       xlabel = "dose efetiva Ang II (ng/kg/min)",
                       xlabelsize = 12,
                       ylabel = "efeito relativo E_AII",
                       ylabelsize = 12)
    xs_aii = collect(Float32, range(0.0, DOSE_MAX, length = 250))
    obs_hill_sem  = Observable(Point2f[(0, 0)])
    obs_hill_ieca = Observable(Point2f[(0, 0)])
    obs_hill_bra  = Observable(Point2f[(0, 0)])
    obs_ec_sem    = Observable([Point2f(0, 0)])
    obs_ec_bra    = Observable([Point2f(0, 0)])
    obs_pt_aii    = Observable([Point2f(0, 0)])
    lines!(ax_hill_aii, obs_hill_ieca; color = AMARELO, linewidth = 5, linestyle = :dot, label = "IECA (r = 1,0)")
    lines!(ax_hill_aii, obs_hill_sem; color = VERDE, linewidth = 2.8, label = "Sem IECA/BRA (r = 1,0)")
    lines!(ax_hill_aii, obs_hill_bra; color = VERMELHO, linewidth = 2.8, label = "BRA crônico (r calibrado)")
    scatter!(ax_hill_aii, obs_ec_sem; color = VERDE, markersize = 9, label = "EC50")
    scatter!(ax_hill_aii, obs_ec_bra; color = VERMELHO, markersize = 9, label = "EC50' (BRA)")
    scatter!(ax_hill_aii, obs_pt_aii; color = :white, strokecolor = LARANJA, strokewidth = 2.5, markersize = 12, label = "instante atual")
    axislegend(ax_hill_aii; position = :lt, framevisible = false, labelsize = 10, patchsize = (16, 10), rowgap = 1)

    sub_rp = Observable("")
    ax_rp_comp = Axis(gg[2, 1];
                      title = "C · Contribuição na Resistência Periférica R_p(t)",
                      titlesize = 13,
                      subtitle = sub_rp,
                      subtitlesize = 11,
                      xlabel = "tempo desde o início da Ang II (min)",
                      xlabelsize = 12,
                      ylabel = "R_p (mmHg·s/ml)",
                      ylabelsize = 12,
                      xticks = -60:30:180)
    obs_rp_base = Observable(Point2f[(0, 0), (1, 0)])
    obs_rp_ne   = Observable(Point2f[(0, 0), (1, 0)])
    obs_rp_tot  = Observable(Point2f[(0, 0), (1, 0)])
    band!(ax_rp_comp, obs_rp_base, obs_rp_ne; color = comalfa(AZUL, 0.30))
    band!(ax_rp_comp, obs_rp_ne, obs_rp_tot; color = comalfa(LARANJA, 0.35))
    lines!(ax_rp_comp, obs_rp_base; color = GRIS, linestyle = :dash, linewidth = 1.8, label = "R_p basal (choque)")
    lines!(ax_rp_comp, obs_rp_ne; color = AZUL, linewidth = 2.2, label = "+ parcela Noradrenalina")
    lines!(ax_rp_comp, obs_rp_tot; color = LARANJA, linewidth = 2.8, label = "R_p total (+ Ang II)")
    vlines!(ax_rp_comp, TSIM; color = (:white, 0.45))
    axislegend(ax_rp_comp; position = :lt, framevisible = false, labelsize = 10, patchsize = (16, 10), rowgap = 1)

    sub_int = Observable("")
    ax_int_dose = Axis(gg[2, 2];
                       title = "D · Interação Temporal: Doses × Resposta da PAM",
                       titlesize = 13,
                       subtitle = sub_int,
                       subtitlesize = 11,
                       xlabel = "tempo desde o início da Ang II (min)",
                       xlabelsize = 12,
                       ylabel = "Ang II (ng/kg/min) / NE×100",
                       ylabelsize = 12,
                       ylabelcolor = LARANJA,
                       yticklabelcolor = LARANJA,
                       xticks = -60:30:180)
    ax_int_pam = Axis(gg[2, 2];
                      yaxisposition = :right,
                      ylabel = "PAM (mmHg)",
                      ylabelsize = 12,
                      ylabelcolor = VERDE,
                      yticklabelcolor = VERDE,
                      backgroundcolor = :transparent)
    hidexdecorations!(ax_int_pam); hidespines!(ax_int_pam)
    hideydecorations!(ax_int_pam; label = false, ticklabels = false, ticks = false)
    linkxaxes!(ax_int_dose, ax_int_pam)

    vlines!(ax_int_dose, [0.0]; color = (:white, 0.35), linestyle = :dash, linewidth = 1.2)
    hlines!(ax_int_pam, [PAM_ALVO_ATHOS]; color = VERDE, linestyle = :dot, linewidth = 1.2)
    hlines!(ax_int_pam, [PAM_ALVO_SSC]; color = GRIS, linestyle = :dot, linewidth = 1.2)

    obs_int_aii  = Observable(Point2f[(0, 0)])
    obs_int_nora = Observable(Point2f[(NaN, NaN)])
    obs_int_pam  = Observable(Point2f[(0, 0)])
    lines!(ax_int_dose, obs_int_aii; color = LARANJA, linewidth = 2.8, label = "Ang II (ng/kg/min)")
    lines!(ax_int_dose, obs_int_nora; color = AZUL, linewidth = 2.8, label = "Nora × 100 (µg/kg/min)")
    lines!(ax_int_pam, obs_int_pam; color = VERDE, linewidth = 3.2, label = "PAM (mmHg)")
    vlines!(ax_int_dose, TSIM; color = (:white, 0.45))
    axislegend(ax_int_dose; position = :rb, framevisible = false, labelsize = 10, patchsize = (16, 10), rowgap = 1)

    rowgap!(gg, 14)
    colgap!(gg, 18)

    function atualizar_curvas_pd!()
        d = DADOS[]; k = K[]
        EmaxN = d.trocas[k].EmaxN
        ecN   = PARAM.ne_ec50[]

        obs_hill_ne[] = Point2f.(xs_ne, Float32.(hill.(xs_ne, EmaxN, ecN)))
        obs_ec50_ne[] = [Point2f(ecN, hill(ecN, EmaxN, ecN))]
        limits!(ax_hill_ne, 0, DOSE_NE_MAX, 0, max(0.35, EmaxN * 1.15))

        y_sem = Float32.(hill.(xs_aii, d.Emax, d.EC50))
        y_bra = Float32.(hill.(xs_aii, d.Emax, d.rb * d.EC50))
        obs_hill_sem[]  = Point2f.(xs_aii, y_sem)
        obs_hill_ieca[] = Point2f.(xs_aii, y_sem)
        obs_hill_bra[]  = Point2f.(xs_aii, y_bra)
        obs_ec_sem[]    = d.EC50 <= DOSE_MAX ? [Point2f(d.EC50, hill(d.EC50, d.Emax, d.EC50))] : [Point2f(NaN, NaN)]
        ec_bra_val      = d.rb * d.EC50
        obs_ec_bra[]    = ec_bra_val <= DOSE_MAX ? [Point2f(ec_bra_val, hill(ec_bra_val, d.Emax, ec_bra_val))] : [Point2f(NaN, NaN)]
        limits!(ax_hill_aii, 0, DOSE_MAX, 0, d.Emax * 1.15)
        return
    end

    function atualizar_series_tratamento!()
        d = DADOS[]; f = FONTE[]; k = K[]
        s = d.series[f][k]
        ii = 1:PASSO_GRAFICO:length(s.t)
        ts = Float32.(s.t[ii])
        x0 = f == :troca ? -FASE_NE : 0.0
        rp0 = Float32(d.plantas[k].Rp)

        e_ne_vec  = s.e_ne === nothing ? zeros(Float32, length(ts)) : Float32.(s.e_ne[ii])
        e_aii_vec = Float32.(s.e_aii[ii])
        rp_base_v = fill(rp0, length(ts))
        rp_ne_v   = rp0 .* (1.0f0 .+ e_ne_vec)
        rp_tot_v  = rp0 .* (1.0f0 .+ e_ne_vec .+ e_aii_vec)

        obs_rp_base[] = Point2f.(ts, rp_base_v)
        obs_rp_ne[]   = Point2f.(ts, rp_ne_v)
        obs_rp_tot[]  = Point2f.(ts, rp_tot_v)

        ymax_rp = max(0.65, Float64(maximum(rp_tot_v)) * 1.25)
        limits!(ax_rp_comp, x0, HORIZONTE, 0.0, ymax_rp)
        sub_rp[] = f == :troca ? "R_p(t) = R_p,0 · (1 + E_NE + E_AII) [modo desmame]" :
                                 "R_p(t) = R_p,0 · (1 + E_AII) [nora fixa já em R_p,0]"

        obs_int_aii[]  = Point2f.(ts, Float32.(s.aii[ii]))
        obs_int_nora[] = s.nora === nothing ? [Point2f(NaN, NaN)] : Point2f.(ts, Float32.(100.0 .* s.nora[ii]))
        obs_int_pam[]  = Point2f.(ts, Float32.(s.pam[ii]))

        ymin_pam = min(55.0, floor(minimum(s.pam) - 2))
        ymax_pam = max(80.0, ceil(maximum(s.pam) + 2))
        limits!(ax_int_dose, x0, HORIZONTE, 0, 88)
        limits!(ax_int_pam, x0, HORIZONTE, ymin_pam, ymax_pam)
        sub_int[] = f == :troca ? "Laranja: Ang II | Azul: Nora×100 | Verde: PAM (eixo dir.)" :
                                  "Laranja: Ang II | Verde: PAM (eixo dir.) — nora fixa"
        return
    end

    function atualizar_pontos_instante!()
        inst = INST[]
        inst === nothing && return
        obs_pt_aii[] = [Point2f(inst.x_aii, inst.e_aii)]
        if inst.x_ne !== nothing && inst.e_ne !== nothing
            obs_pt_ne[] = [Point2f(inst.x_ne, inst.e_ne)]
        else
            obs_pt_ne[] = [Point2f(NaN, NaN)]
        end
        return
    end

    on(_ -> (atualizar_curvas_pd!(); atualizar_series_tratamento!(); atualizar_pontos_instante!()), DADOS)
    on(_ -> (atualizar_series_tratamento!(); atualizar_pontos_instante!()), FONTE)
    on(_ -> (atualizar_curvas_pd!(); atualizar_series_tratamento!(); atualizar_pontos_instante!()), K)
    on(_ -> atualizar_pontos_instante!(), INST)

    atualizar_curvas_pd!()
    atualizar_series_tratamento!()
    atualizar_pontos_instante!()
    colsize!(corpo, 1, Fixed(350))
    return
end