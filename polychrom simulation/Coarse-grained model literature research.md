# Coarse-grained model literature research

## Bottom–up modeling of chromatin segregation due to epigenetic modifications

Lorenzo Boninsegna, Asli Yildirim, Yuxiang Zhan, Frank Alber, Integrative approaches in genome structure analysis, Structure, **30**, 1, (24-36), (2022).https://doi.org/10.1016/j.str.2021.12.003

Model detail：

<img src="C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20221204200114103.png" alt="image-20221204200114103" style="zoom: 50%;" />

- ∼400,000 nucleosomes (of human chromosome 16) as an individual bead.

- choose a confinement 1.8 μm (roughly the size of a chromosome territory)

- interaction

  - electrostatic interactions,
  - van der Waals forces,
  - steric interactions,
  - interactions with a complex solvent containing a zoo of interacting proteins.

- take a coarse-grained approach and define the interaction free energy from the sum of these effects
  $$
  F_{\text {int }}\left(\phi_c\right)= \begin{cases}\chi \Delta^3 \phi_c^2 & \text { if } \phi_{\mathrm{c}}<0.5 \\ \infty & \text { otherwise }\end{cases}
  $$
  

- Monte Carlo Algorithm. Do operation as follows:

  - rotate a segment of polymer about the axis which runs through its ends, 
  - change the binding state of the histone tails,
  - rotate a single bead, translate beads, and pivot the end of the chain. 

  Moves were performed repeatedly to <u>**bring the system to equilibrium**</u>



## Deciphering Phase-Sepertion Mechanism of Double Strnded DNA Induced DNA and Pritein Co-Condensation

Qi Zhi, on Quantitative Biology Symposium 2022

![image-20221204211453268](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20221204211453268.png)