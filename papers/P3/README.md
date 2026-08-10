# A Counterexample to Fourier Alignment in Single-Neuron Modular Addition

**Gautam Neelakantan Memana** · August 2026 · MAIS-P3 · released

**[arXiv:2608.04451](https://arxiv.org/abs/2608.04451)** · [PDF](https://arxiv.org/pdf/2608.04451) · [Provenance](PROVENANCE.md) · resolves [MAIS-O60](../../open-problems/MAIS-O60.md)

## Abstract

We give a negative solution to MAIS-O60. We first construct an example in which an initially active ReLU neuron becomes completely inactive in finite time and thereafter remains frozen at a limit whose Fourier energy is equally distributed among all nonzero real frequency classes. The counterexample holds on an open set of initial conditions and therefore occurs with positive probability under Gaussian initialization. An appendix prepared by GPT-5.6 Sol strengthens the counterexample by showing that the same failure can occur for every Clarke trajectory from an open set of initial conditions, under the convention $\mathrm{ReLU}'(0)=0$, for smooth dead-zone approximations of ReLU, and for fixed-step full-batch gradient descent. Thus, single-frequency alignment is not a general consequence of training a single neuron on modular addition.
