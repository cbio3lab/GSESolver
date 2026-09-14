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

Click here to access our Google Colab Notebook: [https://colab.research.google.com/drive/1bzTD7NvdyLNfh4UF2GJ0UN-swyYTDD4P?usp=sharing](https://colab.research.google.com/drive/1bzTD7NvdyLNfh4UF2GJ0UN-swyYTDD4P?usp=sharing) Otherwise, you can download the notebook from this repo and run it locally. 

---

This repository contains every script, dataset, and supplementary information of our manuscript "Data-driven critical evaluation of the General Solubility Equation", *J. Chem. Inf. Model.*, **2026**, *66 (17)*: 11330–11347.  DOI: https://doi.org/10.1021/acs.jcim.6c01182

<img src="https://acs.silverchair-cdn.com/acs/content_public/journal/jcisd8/issue/66/17/2/jcisd8.2026.66.issue-17.xlargecover-4.jpeg?Expires=1792416707&Signature=5KVZy4whht5LJ3XviKEP5U7TuYcnhk0I5ovbq80TTVAXEudiPMnOVUTNzesdOG3C5AgK-2KXLDE-UpRxf2-aHCl1iqiih4~2drZg8KAHm0CpGnyIbT4FbFk97vexcRXIGVddQsDNTuerrsCHkqHQ6NcUqF0etI4smKnxvHccCSyTHoMDCCccMNqFVHybwp1rh1e0vzCoBMSaaUQMZccqp-aGLGPbZT8WUIOVciryaIar8qk~1t2usAXRrG3nvTUi-FTZvSrrRFprjJh0Bc2tjgo301iVFyL912ZK3TBb85jCOJgZyDdnYln68W~SYwR1aRcCh9XdIGsgG7BBtYN-tg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" alt="A descriptive summary" width="500">


---

**RECENT UPDATES IN OUR COLAB NOTEBOOK:**

* 23/02/2026 *(Added GSESolver logo and description)*
* 01/04/2026 *(Fixed a small RDKit incompatibility)*
* 14/09/2026 *(Added new citation link)*
