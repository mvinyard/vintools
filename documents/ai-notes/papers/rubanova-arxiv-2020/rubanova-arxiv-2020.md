### **Latent ODEs for Irregularly-Sampled Time Series**

***Authors***: *Yulia Rubanova, Ricky T. Q. Chen, David Duvenaud; University of Toronto and the Vector Institute*

### Introduction
- RNNs are bad at irregularly-sampled time series data. A standard trick to circumvent this problem is to divide the irregular time series into equally-sized intervals; one can then impute / aggregate observations using averages. **this approach destroys information, particularly timing, which can be informative of latent variables.**
    - ***note***: *this could be especially relevant argument for methods that seek to gain an understanding of pseudotime.*
- **better approach**: define a **continuous** time model; latent state defined at *all times*.
    - recent progress towards this approach: RNNs w/ continuous dynamics using an exponential decay between observations

- ***Here***: present their new model, **ODE-RNN** to generalize state transitions in RNNs to continuous time dynamics specified by a NN (as is done in *Neural ODEs*).
- demonstrate their model in two *distinct* continuous-time models
    1. standalone autoregressive model
    2. refine latent ODE model from <a href=''>*Chen et al*</a>, using ODE-RNN as a recognition network

***definition***, **Latent ODEs**: define a generative process over time series based on the deterministic evolution of an initial latent state - can be trained as a variational auto-encoder.
