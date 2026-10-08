# symbolic-sensitivity-analysis

Análise de sensibilidade pelo método diferencial com cálculo simbólico (SymPy): primeira e segunda ordem,
funções vetoriais, validação nos vértices e comparação com outros sistemas de álgebra computacional e com
diferenciação automática.

Código usado para obter os resultados do artigo **"Aplicação da Álgebra Computacional na Análise de
Sensibilidade de Modelos Matemáticos"**. Linguagem: Python 3 (testado com 3.12 e SymPy 1.14).

## Instalação

Obrigatório (uma vez):

```
pip install sympy
```

Opcional, só para `medicao_desempenho.py` (o script pula o que não estiver instalado):

```
pip install symengine numpy jax autograd
```

O Maxima (<https://maxima.sourceforge.io>) também é opcional; precisa do comando `maxima` no PATH.

## Arquivos e tabelas do artigo que reproduzem

| Arquivo(s)                           | Execute                        | Gera / reproduz                                      |
|--------------------------------------|--------------------------------|------------------------------------------------------|
| `SSA.py` + `Relatorio.py`            | `python SSA.py`                | `ssa.html`: Tabela 4 (Modelo Exemplo, 1ª ordem)      |
| `SSA2.py` + `Relatorio2.py`          | `python SSA2.py`               | `ssa2.html`: Tabela 14; Me₂ e Mre₂ de P2 e W         |
| `SSA2_FAO56.py` (+ `Relatorio2.py`)  | `python SSA2_FAO56.py`         | `ssa2_fao56.html`: Tabela 10 (FAO-56, n = 8)         |
| `SSA2_parede.py` (+ `Relatorio2.py`) | `python SSA2_parede.py`        | `ssa2_parede.html`: Tabela 11 (parede, n = 10)       |
| `validacao_vertices.py`              | `python validacao_vertices.py` | Tabelas 5, 6, 12 e 13 e as comparações da seção 4.2  |
| `medicao_desempenho.py`              | `python medicao_desempenho.py` | Tabelas 7, 8 e 9 (precisão e tempos dos métodos)     |

Mantenha os arquivos na mesma pasta: cada `SSA*.py` importa o `Relatorio*.py` correspondente.
Os relatórios saem em HTML (vírgula decimal); abra no navegador ou importe no Excel (UTF-8).

## Resultados esperados (conferência rápida)

- `ssa.html`: soma da linha Me = 1,19021; Pie = 0,10257; Mre = 0,58020.
- `ssa2.html`: P2 = 2,05137 e W = 16,81725; Me₂ de P2 = 1,46413 e de W = 4,86863.
- `ssa2_fao56.html`: ET0 = 5,46314; Me total = 0,53889; Mre total = 0,09864.
- `ssa2_parede.html`: q = 18,82564; Me total = 3,35391; Mre total = 0,17816.

## Para configurar outro modelo

Em `SSA.py` ou `SSA2.py`, altere somente: `n`, os valores de `r` e `Δp`, a lista `var`, a linha que
declara os símbolos e a função `f` (ou `F`). As demais linhas não mudam.

## Observações

- Os valores nominais dos modelos FAO-56 e parede são **ilustrativos**, não medições.
- `medicao_desempenho.py` usa, por padrão, os modelos M1, M2, M3, W14, W24 e W54 (cerca de 1 a 2 minutos).
  Para reproduzir as oito colunas das Tabelas 7 a 9:

  ```
  python medicao_desempenho.py --modelos M1,M2,M3,W14,W24,W54,W104,W204 --max-subs 104
  ```

  Leva alguns minutos; o `subs` do SymPy não é executado acima de `--max-subs`, pois o tempo cresce muito.
- Os tempos dependem do computador e das versões (em outro computador chegaram a ser 3 a 5 vezes maiores);
  a ordenação dos métodos se mantém.
- Os desvios do SymEngine podem diferir levemente entre plataformas (arredondamento de ponto flutuante);
  continuam na ordem de 1e-15 ou menos.
- Diferenças finitas: a linha "Dif. finitas" usa escalares NumPy (como no artigo); "Dif. finitas (float)"
  usa ponto flutuante nativo e é de 3 a 8 vezes mais rápida. O tempo desse método depende da implementação.
- No Windows, se `python` não funcionar, use `py`. Salve os `.py` em UTF-8 se for editá-los.

## Como citar

Se usar este código, cite o artigo "Aplicação da Álgebra Computacional na Análise de Sensibilidade de
Modelos Matemáticos" (referência completa a ser acrescentada após a publicação).