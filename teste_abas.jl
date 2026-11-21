# Protótipo das abas: só testa a navegação (ainda sem o modelo).
# Cinco abas, com um estado compartilhado (slider, relógio e pausa) que NÃO pode reiniciar ao trocar de aba.
#
#   julia teste_abas.jl tela      cada aba é uma figura; a mesma janela troca de figura
#   julia teste_abas.jl ocultar   uma figura só; as abas que não estão ativas ficam ocultas
using GLMakie

set_theme!(theme_dark())

const NOMES   = ["Monitoramento", "Tratamento", "Modelo 3D", "Cenários", "Engenharia"]
const FUNDO   = RGBf(0.07, 0.09, 0.13)
const ATIVA   = RGBf(0.38, 0.64, 1.00)
const INATIVA = RGBf(0.20, 0.25, 0.34)

# estado compartilhado entre as abas
const K       = Observable(5.0)       # valor do slider
const T       = Observable(0.0)       # relógio da "simulação" (min)
const TOCANDO = Observable(true)
const ULTIMO  = Ref(time())

# o relógio anda com o tempo real, mesmo se a aba ficou escondida por um tempo
function avancar!()
    agora = time()
    TOCANDO[] && (T[] += 2.0 * (agora - ULTIMO[]))
    ULTIMO[] = agora
    return
end

# barra com os cinco botões; devolve a lista de botões
function barra!(gl, ativa, trocar!)
    bts = Any[]
    for j in 1:5
        bt = Button(gl[1, j]; label = NOMES[j], width = 190,
                    buttoncolor = j == ativa ? ATIVA : INATIVA,
                    labelcolor = j == ativa ? :black : :white)
        on(_ -> trocar!(j), bt.clicks)
        push!(bts, bt)
    end
    return bts
end

# conteúdo de uma aba; `reg` recebe os blocos criados (usado para ocultar)
function conteudo!(gl, i, reg)
    r(b) = (push!(reg, b); b)
    ax = r(Axis(gl[1, 1]; title = "Aba $i: $(NOMES[i])", xlabel = "tempo (min)", ylabel = "valor",
                limits = (0, 60, -11, 11)))
    xs = 0:0.2:60
    lines!(ax, lift(k -> Point2f.(xs, k .* sin.(xs ./ (1 + i))), K); linewidth = 3, color = ATIVA)
    scatter!(ax, lift(t -> [Point2f(mod(t, 60), 0)], T); markersize = 14, color = :orange)
    r(Label(gl[2, 1], lift((t, k) -> "Relógio da simulação: $(round(t, digits = 1)) min   |   controle: $(round(k, digits = 1))", T, K)))
    sl = r(Slider(gl[3, 1]; range = 1:0.5:10, startvalue = K[], width = 400))
    on(sl.value) do v
        v != K[] && (K[] = v)
    end
    on(K) do v
        sl.value[] != v && set_close_to!(sl, v)
    end
    bp = r(Button(gl[4, 1]; label = TOCANDO[] ? "Pausar" : "Reproduzir", width = 150))
    on(_ -> (TOCANDO[] = !TOCANDO[]), bp.clicks)
    on(TOCANDO) do p
        bp.label[] = p ? "Pausar" : "Reproduzir"
    end
    return
end

# ---------------------------------------------------------------- modo 1: a janela troca de figura

function modo_tela()
    figs = Vector{Any}(undef, 5)
    tela = Ref{Any}(nothing)
    trocar!(j) = display(tela[], figs[j])
    for i in 1:5
        fig = Figure(size = (1100, 650), backgroundcolor = FUNDO)
        gb = fig[1, 1] = GridLayout()
        barra!(gb, i, trocar!)
        gc = fig[2, 1] = GridLayout()
        conteudo!(gc, i, Any[])
        on(events(fig).tick) do _
            avancar!()
        end
        figs[i] = fig
    end
    tela[] = display(figs[1])
    isinteractive() || wait(tela[])
    return
end

# ---------------------------------------------------------------- modo 2: uma figura, abas ocultas

function visivel!(b, v)
    b.blockscene.visible[] = v
    try
        b.scene.visible[] = v
    catch
    end
    return
end

function modo_ocultar()
    fig  = Figure(size = (1100, 650), backgroundcolor = FUNDO)
    regs = [Any[] for _ in 1:5]
    bts  = Any[]
    function mostrar!(j)
        for i in 1:5, b in regs[i]
            visivel!(b, i == j)
        end
        for (n, bt) in enumerate(bts)
            bt.buttoncolor[] = n == j ? ATIVA : INATIVA
            bt.labelcolor[]  = n == j ? :black : :white
        end
        return
    end
    gb = fig[1, 1] = GridLayout()
    append!(bts, barra!(gb, 1, mostrar!))
    for i in 1:5
        gc = fig[2, 1] = GridLayout()
        conteudo!(gc, i, regs[i])
    end
    mostrar!(1)
    on(events(fig).tick) do _
        avancar!()
    end
    tela = display(fig)
    isinteractive() || wait(tela)
    return
end

modo = isempty(ARGS) ? "tela" : lowercase(ARGS[1])
if modo == "tela"
    modo_tela()
elseif modo == "ocultar"
    modo_ocultar()
else
    error("Use: julia teste_abas.jl tela   ou   julia teste_abas.jl ocultar")
end