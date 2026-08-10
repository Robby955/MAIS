# Progress on MAIS-O62

*Progress report for [MAIS-O62](../MAIS-O62.md): Minimal network width for exact modular addition · Opus and Iris Shi · August 2026 · Status: verified; MAIS-O62 remains open.*

## Current draft

**[Read the current draft (PDF)](https://github.com/user-attachments/files/30774709/MAIS_O62_Bound_Version_2.pdf)**

Revision history, supporting files, and discussion are available in [issue #2](https://github.com/lionellevine/MAIS/issues/2).

## Width lower bound and the case $p=3$

Opus and Iris Shi prove the improved lower bound

$$
H_{\min}(p)\ge \frac{3p-1}{2}
$$

for every odd integer $p$, and give an exact width-4 construction establishing $H_{\min}(3)=4$.

## Status of the remaining claims

- At $p=5$, the supplied width-8 fits are numerical candidates with residuals near machine precision. They do not prove $H_{\min}(5)\le8$.
- At $p=7$, the Fourier construction gives an exact width-13 fit. The separate numerical width-13 candidate adds no upper bound.
- Searches found no width-7 fit at $p=5$ and no width-12 fit at $p=7$. These negative searches are documented and calibrated, but do not prove nonexistence.

*Related: [MAIS-O62](../MAIS-O62.md) (the open problem) · [issue #2](https://github.com/lionellevine/MAIS/issues/2) (submission and discussion).*
