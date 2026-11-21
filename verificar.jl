# Confere o modelo contra os números reportados no documento:  julia verificar.jl
include(joinpath(@__DIR__, "modelo.jl"))
using Printf

const GANHO = 1.0   # mesmo valor inicial do aplicativo

println("== 1. Plantas (documento: tabela 'Entradas, parâmetros e saídas') ==")
@printf("%-16s %5s %6s %6s %6s %6s %7s | %-10s %4s %4s | %s\n",
        "Cenário", "FC", "VS", "DC", "T_ej", "R_p", "Z_c", "PAS/PAD", "PAM", "PP", "documento")
doc = ["124/74, PAM 95, PP 50", "145/68, PAM 102, PP 77", "81/38, PAM 58, PP 44", "103/29, PAM 62, PP 74"]
for (p, d) in zip((planta_adulto(), planta_idoso(), planta_choque_adulto(), planta_choque_idoso()), doc)
    b = batimento(p)
    @printf("%-16s %5.0f %6.1f %6.2f %6.3f %6.3f %7.4f | %4.0f/%-5.0f %4.0f %4.0f | %s\n",
            p.nome, p.FC, p.VS, dc_lmin(p), t_ejecao(p.FC), p.Rp, zc(p.Rp), b.PAS, b.PAD, b.PAM, b.PP, d)
end
a = planta_adulto()
@printf("Checagem cruzada: C ajustado %.2f contra VS/PP = %.0f/50 = %.2f ml/mmHg (documento: 1,91 e 1,88)\n",
        a.C, a.VS, a.VS / 50)

println("\n== 2. Calibração da Ang II (documento: r = 11,9; faixa 10,5 a 14,8) ==")
for Emax in (0.3, 0.5, 1.0)
    @printf("Emax %.1f: EC50 = %6.1f ng/kg/min, r do BRA = %.1f\n", Emax, ec50_aii(Emax), r_bra(Emax))
end

println("\n== 3. Comparativo em 3 h (ganho do passo = $GANHO ng/kg/min por mmHg) ==")
@printf("%-24s %5s %5s | %6s %6s %6s %6s %9s %8s %6s\n", "Grupo", "Emax", "r",
        "média", "final", "PAM1h", "PAM3h", "PAS/PAD", "atinge", "F15")
grupos = [("Adulto", planta_choque_adulto(), false), ("Idoso sem IECA/BRA", planta_choque_idoso(), false),
          ("Idoso com BRA crônico", planta_choque_idoso(), true)]
for (nome, p, bra) in grupos, Emax in (0.5, 0.3, 1.0)
    r = bra ? r_bra(Emax) : 1.0
    s = titular(p, Emax, r; ganho = GANHO)
    b = batimento(com_rp(p, s.rp_final))
    alvo = s.t_alvo === nothing ? "não" : @sprintf("%.0f min", s.t_alvo)
    @printf("%-24s %5.1f %5.1f | %6.1f %6.1f %6.1f %6.1f %4.0f/%-4.0f %8s %6.2f\n",
            nome, Emax, r, s.dose_media, s.dose_final, s.pam_1h, s.pam_3h, b.PAS, b.PAD, alvo, s.F15)
end
println("Documento (Emax 0,5; faixa 0,3 a 1,0):")
println("  Adulto:   média 65 (50 a 72), final 77 (54 a 80), PAM 1 h 74 (71 a 75), 3 h 75 (100/53)")
println("  Idoso:    média 36 (33 a 46), final 39 (34 a 55), PAM 1 h 75 (74 a 75), 3 h 75 (116/40), entre 36 e 85 min")
println("  Idoso BRA: média 76, final 80 saturada, PAM 65,5 (106/32), não atinge; F15 ~0,77 sem uso e ~0,16 com BRA")

println("\n== 4. Aproximação quase-estática contra Euler 0,1 ms nas 3 h (idoso sem IECA/BRA, Emax 0,5) ==")
function euler_completo(p::Planta, Emax, r; ganho)
    EC50 = ec50_aii(Emax);  k = log(2) / (T_MEIA_AII * 60);  dt = DT_EULER
    T = 60 / p.FC;  Tej = t_ejecao(p.FC)
    Pc = p.Rp * dc_mls(p);  u = DOSE_INICIAL;  x = 0.0
    tc = 0.0;  soma = 0.0;  nb = 0;  pam_bat = pam_regime(p)
    prox = INTERVALO_AJUSTE * 60;  pam1h = NaN
    for i in 1:round(Int, HORIZONTE * 60 / dt)
        t  = i * dt
        Rp = p.Rp * (1 + hill(x, Emax, r * EC50))
        q  = q_in(tc, p.VS, Tej)
        soma += Pc + zc(Rp) * q;  nb += 1
        Pc += dt * (q - (Pc - P_VEN) / Rp) / p.C
        x  += dt * k * (u - x)
        tc += dt
        if tc >= T
            pam_bat = soma / nb;  soma = 0.0;  nb = 0;  tc -= T
        end
        if t >= prox
            u = clamp(u + clamp(ganho * (PAM_ALVO_ATHOS - pam_bat), -PASSO_MAX, PASSO_MAX), 0.0, DOSE_MAX)
            prox += INTERVALO_AJUSTE * 60
        end
        isnan(pam1h) && t >= 3600 && (pam1h = pam_bat)
    end
    return (pam_1h = pam1h, pam_3h = pam_bat, dose_final = u)
end
p = planta_choque_idoso()
e = euler_completo(p, 0.5, 1.0; ganho = GANHO)
s = titular(p, 0.5, 1.0; ganho = GANHO)
@printf("Euler 0,1 ms:   PAM 1 h %.2f, PAM 3 h %.2f, dose final %.2f\n", e.pam_1h, e.pam_3h, e.dose_final)
@printf("Quase-estática: PAM 1 h %.2f, PAM 3 h %.2f, dose final %.2f\n", s.pam_1h, s.pam_3h, s.dose_final)
