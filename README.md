# symbolic-sensitivity-analysis

Sensitivity analysis using the differential method with symbolic computation (SymPy): first and second order, vector functions, vertex validation, and comparison with other computer algebra systems and automatic differentiation.

Code used to obtain the results for the article **"Application of Computer Algebra in the Sensitivity Analysis of Mathematical Models"**. Language: Python 3 (tested with 3.12 and SymPy 1.14).

## Installation

Required (once):

```
pip install sympy

```

Optional, only for `medicao_desempenho.py` (the script skips what is not installed):

```
pip install symengine numpy jax autograd

```

Maxima ([https://maxima.sourceforge.io](https://maxima.sourceforge.io)) is also optional; it requires the `maxima` command in your PATH.

## Files and the article tables they reproduce

| File(s) | Run | Generates / reproduces |
| --- | --- | --- |
| `SSA.py` + `Relatorio.py` | `python SSA.py` | `ssa.html`: Table 4 (Example Model, 1st order) |
| `SSA2.py` + `Relatorio2.py` | `python SSA2.py` | `ssa2.html`: Table 14; Me₂ and Mre₂ of P2 and W |
| `SSA2_FAO56.py` (+ `Relatorio2.py`) | `python SSA2_FAO56.py` | `ssa2_fao56.html`: Table 10 (FAO-56, n = 8) |
| `SSA2_parede.py` (+ `Relatorio2.py`) | `python SSA2_parede.py` | `ssa2_parede.html`: Table 11 (wall, n = 10) |
| `validacao_vertices.py` | `python validacao_vertices.py` | Tables 5, 6, 12, and 13 and comparisons from Sec. 4.2 |
| `medicao_desempenho.py` | `python medicao_desempenho.py` | Tables 7, 8, and 9 (methods precision and execution times) |

Keep the files in the same folder: each `SSA*.py` imports its corresponding `Relatorio*.py`. The reports are generated in HTML (using a decimal comma); open them in a browser or import them into Excel (UTF-8).

## Expected results (quick check)

* `ssa.html`: sum of the Me row = 1.19021; Pie = 0.10257; Mre = 0.58020.
* `ssa2.html`: P2 = 2.05137 and W = 16.81725; Me₂ of P2 = 1.46413 and of W = 4.86863.
* `ssa2_fao56.html`: ET0 = 5.46314; Total Me = 0.53889; Total Mre = 0.09864.
* `ssa2_parede.html`: q = 18.82564; Total Me = 3.35391; Total Mre = 0.17816.

## Configuring a different model

In `SSA.py` or `SSA2.py`, change only: `n`, the values of `r` and `Δp`, the `var` list, the line declaring the symbols, and the `f` (or `F`) function. The remaining lines do not change.

## Notes

* The nominal values for the FAO-56 and wall models are **illustrative**, not actual measurements.
* `medicao_desempenho.py` uses models M1, M2, M3, W14, W24, and W54 by default (takes about 1 to 2 minutes).
To reproduce the eight columns of Tables 7 to 9:
```
python medicao_desempenho.py --modelos M1,M2,M3,W14,W24,W54,W104,W204 --max-subs 104

```


This takes a few minutes; SymPy's `subs` is not executed above `--max-subs`, as execution time grows significantly.
* Execution times depend on your computer and software versions (on another machine they were up to 3 to 5 times slower); however, the performance ranking of the methods remains consistent.
* SymEngine deviations may differ slightly between platforms (due to floating-point rounding); they remain on the order of 1e-15 or less.
* Finite differences: the "Dif. finitas" row uses NumPy scalars (as in the article); "Dif. finitas (float)" uses native floating-point operations and is 3 to 8 times faster. The speed of this method depends heavily on the implementation.
* On Windows, if the `python` command doesn't work, try using `py`. Save your `.py` files in UTF-8 if you edit them.

## How to cite

If you use this code, please cite the article "Application of Computer Algebra in the Sensitivity Analysis of Mathematical Models" (full reference to be added after publication).ão).
