# Modelo_cardiovascular_PBL_Felipe
Mathematical modeling and simulation of refractory vasodilatory shock in elderly patients using a three-element Windkessel model and vasoactive drugs.

## Como rodar

Requer Julia 1.10 ou mais recente. Na pasta do projeto, instale as dependências uma vez:

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Depois:

```
julia --project=. verificar.jl     # confere o modelo contra os números do documento (texto no terminal)
julia --project=. app_abas.jl      # abre a plataforma em abas (GLMakie)
```

## Arquivos

- `modelo.jl`: modelo Windkessel de 3 elementos + PK/PD da angiotensina II e da noradrenalina
- `verificar.jl`: comparação dos resultados com o documento
- `app_abas.jl`: aplicativo principal (abas 1 e 4) e estado compartilhado
- `aba_tratamento.jl`, `aba_modelo3d.jl`, `aba_engenharia.jl`: abas 2, 3 e 5
- `teste_abas.jl`: protótipo antigo da navegação entre abas (não usado pelo app)
