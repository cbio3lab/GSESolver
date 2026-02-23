# **GSESolver**
A FAIR (*Findable, Accessible, Interoperable, and Reusable*) workflow to the accuracy of the General Solubility Equation for small molecules.
---

<div align="center;">
<img src="https://i.imgur.com/BSUS5ck.png" width="800">
</div>

---

Intrinsic solubility (expressed as $\log{S_0}$) indicates the maximum concentration of a compound in its neutral state that can be dissolved in water. This physicochemical property is evergrowingly assessed in drug design/developement, medicinal chemistry, and agrochemical sciences. 

It is often not possible to experimentally measure the $\log{S_0}$. Thus, the General Solubility Equation (**GSE**) is oftentimes used for a rapid assessment, which uses the molecule's melting point ($mp$ in °C) and *n*-octanol/water partition coefficient ($\log{P_{\text{N}}}$):

$$\log{S_0}=0.5-0.01(mp-25)-\log{P_{\text{N}}}$$

*When should I use the GSE?* This notebook provides an Artificial Neural Network model that assesses if your molecules of interest might have an accurate $\log{S_0}$ prediction using the GSE. 

Click here to access our Google Colab Script: [https://colab.research.google.com/drive/1bzTD7NvdyLNfh4UF2GJ0UN-swyYTDD4P?usp=sharing](https://colab.research.google.com/drive/1bzTD7NvdyLNfh4UF2GJ0UN-swyYTDD4P?usp=sharing)

---

This repository contains every script, dataset, and supplementary information of our manuscript "Critical evaluation of the General Solubility Equation: When is it a valid solubility predictor?", *Working paper*, **2026**.

---

**RECENT UPDATES IN OUR COLAB NOTEBOOK:**

* 23/02/2026 *(Added GSESolver logo and description)*
