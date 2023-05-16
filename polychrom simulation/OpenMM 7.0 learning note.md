# Learning Note of OpenMM 7.0

## previous work using OpenMM dynamics simulation toolkit

Rajpurkar, A.R. et al., *Nature Communication* , 2021

 https://doi.org/10.1038/s41467-021-23831-4

![image-20220522170859370](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220522170859370.png)

首先，就是之前调研的这篇深度学习的工作。（以距离矩阵为输入，转录开关为学习对象。）其中，使用了OpenMM来模拟一些简单的带enhancer & promoter的polymer（较为简单的转录开关机制），来判断一些判断方法的实用性。模拟中，使用了下文中对于OpenMM的改进，即加入了force calculation functions，使得本文模拟中，52个结合位点间没有吸引，只是靠着density提供吸引力（用来模拟核约束），而在位点过于接近时候有一个斥力。

Nuebler, J. et al.,*Proceedings of the National Academy of Sciences* , 2018

https://www.pnas.org/doi/suppl/10.1073/pnas.1717730115

![image-20220522170643875](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220522170643875.png)

这篇是上一篇的基础，也使用这一系统。本文主要研究哺乳动物染色质上TADs，并使用分子链模拟这一结构的形成，以及一些生物分子浓度改变时，TADs会发生何种变化。

<img src="C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220522204823576.png" alt="image-20220522204823576" style="zoom: 50%;" />

[Supplemental Information](C:\Users\zhouquan\Zotero\storage\LZG4WFCS\pnas.1717730115.sapp.pdf)中详细描述了模拟过程。聚合物表示为具有一系列单体，用谐波势能连接，有体积排斥势，并且在两个typeB的单体间的微小吸引力

## basic functions of OpenMM 7.0

- 读取输入文件（Amber, CHARMM, Gromacs, and Desmond）
- 编辑分子模型
- 定义作用在分子系统上的力（明确指定，或从文件加载力场）
- 计算力和能量
- 运行模拟
- 输出结果

其中，我们最关心的肯定是自定义力场部分，这里OpenMM提供的有：

![image-20220522185054095](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220522185054095.png)

## simulation test

首先测试[教程](http://docs.openmm.org/latest/userguide/application/02_running_sims.html#a-first-example)自带的[模拟程序](C:\Users\zhouquan\anaconda3\Library\share\openmm\examples)。通过这个程序我们可以对模拟环境进行调整。

首先是开始import的功能包：

```python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
```

运行需要`cd`到对应母文件夹，使用`python simulatePdb.py`命令启动。

### A pdb example

之后导入已经写好分子信息，这里使用的是pdb文件，实际上还有其他的种类，如(Amber和Gromacs)，包含了分子拓扑结构和原子位置

```python
pdb = PDBFile('input.pdb')
```

之后就是引入力场，力场信息使用.xml文件存储，这里引入两个

```python
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
```

然后就是设置模拟系统，主要是加入分子信息和力场，其中有三个关键参数：远距离静电作用为particle mesh Ewald，直接空间作用（应该是碰撞体积）为1nm，约束键长度为氢分子的键长。

```python
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
```

下一步调整的是积分器，这里的示例是郎之万中值积分器（不同的模型应当选取不同的积分器），同样三个参数：温度，摩擦因数和<u>**步长**</u>：

```python
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
```

之后就是启动模拟了，需要三个部分：分子信息，系统信息和积分器，之后就是设置初始位置，然后设置形成一个附近能量极小值，一般来说开始模拟的时候要加上这一步：

```python
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
```

最后就是模拟的输出了，首先会输出一个叫output的文件，每1000步写入一次信息。同时，而有一个在cmd（指定stdout作为输出文件即可）中print的内容（=True的内容会被打印），最后指定模拟步数。

```python
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000,step=True,potentialEnergy=True, temperature=True))
simulation.step(10000)
```

最终效果为：

![image-20220526094659471](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220526094659471.png)

### Using AMBER Files

Amber系统需要导入两个文件，prmtop文件包括了分子拓扑形式和力场参数，inpcrd是原子位置：

```python
prmtop = AmberPrmtopFile('input.prmtop')
inpcrd = AmberInpcrdFile('input.inpcrd')
```

同时，在设置部分，还要额外从inpcrd设置Box Vector

```python
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
```

### Using Gromacs Files

与上部分基本相同，只是replacing `AmberInpcrdFile` and `AmberPrmtopFile` with `GromacsGroFile` and `GromacsTopFile`

### Using CHARMM Files

同样，导入信息也有一定区别。

### The OpenMM-Setup Application

实际上，这一运行脚本可以自动生成，通过openMM setup。安装和启动可以使用以下命令：

```python
conda install -c conda-forge openmm-setup
openmm-setup
```

会生成一个网页，从而根据输入文件生成脚本

### 模拟参数 Simulation Parameters

#### 平台

平台使用的应该是OpenCL, CUDA, CPU, or Reference

#### 力场

这里只大部分请参考文档，只记录一小部分：

amber14提供了一系列标准力场

![image-20220529194244232](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220529194244232.png)

James A. Maier, Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling. Ff14sb: improving the accuracy of protein side chain and backbone parameters from ff99sb. *Journal of Chemical Theory and Computation*, 11(8):3696–3713, 2015.

这里大部分力场文件收录在工具包中，可以使用导入

```python
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
```

####  Constraints

有四个值可以使用

![image-20220529195049092](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220529195049092.png)

#### 积分器 Integrators

其中值得注意的是，有一个布朗积分器，可以用于模拟扩散动力学系统

## Polymer Dynamics

### 之前引用文献的模拟模型

Nuebler, J. et al.,*Proceedings of the National Academy of Sciences* , 2018

https://www.pnas.org/doi/suppl/10.1073/pnas.1717730115

![image-20220522170643875](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220522170643875.png)

首先我们来看本文的模型，这里取DNA是一条有谐波势连接的单体链，这里主要分析一下他的具体模拟包——[polychrom包](https://github.com/open2c/polychrom)。

（以下超链接可以跳转到网页详细说明）在这个包中，力场在 [`polychrom.forces`](https://polychrom.readthedocs.io/en/latest/polychrom.forces.html#module-polychrom.forces) 和 [`polychrom.forcekits`](https://polychrom.readthedocs.io/en/latest/polychrom.forcekits.html#module-polychrom.forcekits) 中定义，接触关系由 [`polychrom.contactmaps`](https://polychrom.readthedocs.io/en/latest/polychrom.contactmaps.html#module-polychrom.contactmaps)文件定义，加载和存储则在[`polychrom.polymerutils`](https://polychrom.readthedocs.io/en/latest/polychrom.polymerutils.html#module-polychrom.polymerutils)中，最终由[`polychrom.simulation`](https://polychrom.readthedocs.io/en/latest/polychrom.simulation.html#module-polychrom.simulation)统合，分析则有[`polychrom.polymer_analyses`](https://polychrom.readthedocs.io/en/latest/polychrom.polymer_analyses.html#module-polychrom.polymer_analyses)脚本负责。





### amber14中的DNA模型

Maier JA et al.,*Journal of chemical theory and computation* ,2015

这一分子模型为原子尺度，我们难以使用。

### polychrom

#### 安装

环境安装：

- [x] conda安装
- [ ] CUDA安装

前置包安装：

- [x] six
- [x] cython
- [x] numpy>=1.9
- [x] scipy>=0.16
- [x] h5py>=2.5
- [x] pandas>=0.19
- [x] joblib
- [x] pyknotid
- [x] pytest

#### 测试

详见滴答清单上操作流程

#### 正式模拟

E-P1系列：

polymer参数设置

```python
N=400

platform="CUDA",
integrator="variableLangevin",
error_tol=0.003,#<0.01
GPU="0",#0 represent for the first GPU, and 1 for the second. If only 1 GPU is available, it will be selected automatically
collision_rate=0.03,#0.01-0.05
N=N,
save_decimals=2,
PBCbox=False,
reporters=[reporter],
temperature=293,#devault value is 300K
mass=1518500,#estimated by dsDNA 607.4 per bp and 2.5Kb per monomer

polymer = starting_conformations.grow_cubic(400, 1000)#(m, n)represents a m-long polymer in a cube of n*n*n
```

Force设置

首先，限制区域的力为：spherical_confinement，Constrain particles to be within a sphere. 

在Harmonic Bonds中，两个参数分别为键长(bL bondLength)和键长波动(bWD bondWiggleDistance)（自由能导致，Bond energy equal to 1kT at bondWiggleDistance）

在angle force中，k代表了聚合物的僵硬程度——越小越软，越大越僵硬，推荐值为1.5。（i-th triplet Potential is k * alpha^2 * 0.5 * kT）

最后，排斥力的参数代表了r=0处的力场势能(It has the value of trunc at zero, stays flat until 0.6-0.7 and then drops to zero together with its first derivative at r=1.0.)

当然，根据我们的现象，需要一个吸引力场，形式可以尝试（可能是本次重点）

最终数据的参数如下表，每次进行100次模拟，模拟长度为5000步长

| Number | 步长 | 键力bL/bWD | 角度力k | 排斥力trunc | 吸引力场 |      |
| ------ | ---- | ---------- | ------- | ----------- | -------- | ---- |
| 1      | auto | 1.0/0.05   | 1.5     | 3.0         | No       |      |
| 2      |      | 0.5/0.025  |         |             |          |      |
| 3      |      |            |         |             |          |      |

步长实际上也是我们讨论的重点，理论上variableLangevin会自定义步长

#### 数据处理：多次实验

但我们发现了一些问题：程序自带的reporter并不能记录所有时刻的monomer位置，仅能汇报最后时刻。所以，目前有两种改进方法，本部分记录第一种。

特别地：<u>**这几次模拟使用的monomers间距为0.5单位**</u>

进行10次持续时间不同的模拟，每一次模拟50个分子链，得到不同时间的数据。按照统计规律来说，此数据得到的MSD随时间分布，应当与实际记录分子链各个时刻数据统计上一致。

进行数据处理后，实验效果不佳。MSD随着时间分布如下：

![image-20220821202552795](C:\Users\zhouquan\AppData\Roaming\Typora\typora-user-images\image-20220821202552795.png)

#### 数据处理：单次模拟

而正常的单次模拟给出的结果如下。我们可以看到，如果不重新建立一个simulation，仅仅`do_block`其不会进行初始化，以下为验证：

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\test_E-P2\位点连续性验证.jpg" alt="位点连续性验证" style="zoom:33%;" />

最终给出的MSD如下图，为500长度分子链中第250个和350个之间相对位移的MSD

![MSD随时间变化](C:\Users\zhouquan\OneDrive\research\polychrom simulation\test_E-P2\MSD随时间变化.jpg)

其中d为两个monomers之间的平均距离。目前有几个问题：

- 时间尺度太小（待使用郎之万积分法中给出的确定）
- 为何z方向上差别如此大？
- 100个monomers之间为何距离如此接近？可能和外加的总体限制立场有关（有待调整）
- 波动太大（这组数据仅为10次实验）

使用VMD三维作图





#### 确定初始状态

##### cubic standard mode

grows a ring starting with a 4-monomer ring in the middle, on a cubic lattice in the cubic box of size boxSize.

| cubic size | figure                                                       |
| ---------- | ------------------------------------------------------------ |
| 500        | <img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\cubic 500.jpg" alt="cubic 500" style="zoom:25%;" /> |
| 1000       | <img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\cubic 1000.jpg" alt="cubic 1000" style="zoom:25%;" /> |
| 2000       | <img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\cubic 2000.jpg" alt="cubic 2000" style="zoom:25%;" /> |

##### cubic linear mode

it grows a linearly organized chain from 0 to size.

| box size |                                                              |
| -------- | ------------------------------------------------------------ |
| 200      | <img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\cubic linear 200.jpg" alt="cubic linear 200" style="zoom:25%;" /> |
| 400      | <img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\cubic linear 400.jpg" alt="cubic linear 400" style="zoom:25%;" /> |

##### spiral

 Creates a "propagating spiral", often used as an easy mitotic-like

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\sprial 5 6.jpg" alt="sprial 5 6" style="zoom:33%;" />

##### random walk

Creates a freely joined chain of length N with step step_size

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\random walk step 1.jpg" alt="random walk step 1" style="zoom: 24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\initial_state_figure\random walk step 1 dot.jpg" alt="random walk step 1 dot" style="zoom: 24%;" />

P.S. constrained_random_walk ?

#### 关于步长

在咨询polychrom作者后，得到了如下建议：

- running a simulation for a given amount of time (in many replicates), calculating average MSD(t)
- fit it to the biologically relevant MSD measured in some experiment
- For all practical purposes it is safe to keep integrator as variableLangevin, and basically measure time in timesteps.

再次询问如何具体将实验与模拟的MSD(t)联系起来，是简单地加两个参数？

### 正式测试

#### 拟合方法

首先，在力场上，分为四项：

- Harmonic Bonds，表示monomers的连接情况
- 角度力angle force代表polymer的僵硬程度
- 排斥力用来模拟碰撞
- 需要一个<u>**吸引力场**</u>，形式可以尝试

其次，在步长确定上，分两种模拟手段：

- variableLangevin分析器，使用timestep或记录时间为t
- Langevin分析器，使用固定timestep（与variableLangevin分析器相差一个量级左右）

在两种手段的模拟结果($MSD_S(t)$)上，将$a\cdot MSD_S(bt)$与$MSD_E(t)$进行对比，从而找到最适合的吸引力场

#### 长度时间比例$a,b$确定

实际上，长度比例比较容易确定：我们使用100个monomers模拟长度为47kbp的一段DNA链，也就是说相邻monomers距离$1A.U.=0.47kbp=0.47\times340nm=160nm$
$$
a=(160nm)^2=2.56\times10^{-14}
$$
而时间这一问题就较为复杂。实际上，理想情况下，不同的步长缩放到同样比例后，应当有近似的结果。如步长为$10fs$（1倍variableLangevin积分器），或$100fs$（10倍variableLangevin积分器），得到相同的结果仅需要改变力场大小（理论上力量纲为$kg\cdot m/s^2$）。理想情况下，放大10倍步长，减小100倍力场大小后情况应该不变。

下面列举步长不同情况下的三个结果（无外力场）



<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\1vLt No appled forcefield.jpg" alt="1vLt No appled forcefield" style="zoom: 24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\10vLt No appled forcefield.jpg" alt="10vLt No appled forcefield" style="zoom: 24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\100vLt No appled forcefield.jpg" alt="100vLt No appled forcefield" style="zoom: 33%;" />



#### 力场确定

##### ~~sphereWell球形吸引势场，$r=20A.U.,depth=1kT$~~

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\10vLt SphereWell forcefield.jpg" alt="10vLt SphereWell forcefield" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\100vLt SphereWell forcefield.jpg" alt="100vLt SphereWell forcefield" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\1000vLt SphereWell forcefield.jpg" alt="1000vLt SphereWell forcefield" style="zoom:24%;" />

*这里参数使用错误，并且势场理解有误，势能为球形，而不是仅仅是势场为球形。



##### 中心势场，中心弹簧势场

1. 首先观察势能系数变化的影响

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.01spring1e3vLt.jpg" alt="0.01spring1e3vLt" style="zoom: 24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.02spring1e3vLt.jpg" alt="0.02spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.03spring1e3vLt.jpg" alt="0.03spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.05spring1e3vLt.jpg" alt="0.05spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.1spring1e3vLt.jpg" alt="0.1spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\1spring1e3vLt.jpg" alt="1spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\noForce1e3vLt.jpg" alt="noForce1e3vLt" style="zoom: 24%;" />

这里最后附上无势场的MSD情况

与实验数据比较（将模拟的xyz平均）

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\compare0.01spring1e3vLt.jpg" alt="compare0.01spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\compare0.02spring1e3vLt.jpg" alt="compare0.02spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\compare0.03spring1e3vLt.jpg" alt="compare0.03spring1e3vLt" style="zoom:25%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\compare0.1spring1e3vLt.jpg" alt="compare0.1spring1e3vLt" style="zoom:24%;" />

2. 其次考虑长度带来的影响

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.01spring1e3vLt.jpg" alt="0.01spring1e3vLt" style="zoom: 24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\long0.01spring1e3vLt.jpg" alt="long0.01spring1e3vLt" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.02spring1e3vLt.jpg" alt="0.02spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\long0.02spring1e3vLt.jpg" alt="long0.02spring1e3vLt" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.1spring1e3vLt.jpg" alt="0.1spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\long0.1spring1e3vLt.jpg" alt="long0.1spring1e3vLt" style="zoom:24%;" />

3. 然后考虑和enhancer和promoter位置的影响

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.01spring1e3vLt.jpg" alt="0.01spring1e3vLt" style="zoom: 24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\expand0.01spring1e3vLt.jpg" alt="expand0.01spring1e3vLt" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\prolong0.01spring1e3vLt.jpg" alt="prolong0.01spring1e3vLt" style="zoom:25%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.02spring1e3vLt.jpg" alt="0.02spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\expand0.02spring1e3vLt.jpg" alt="expand0.02spring1e3vLt" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\prolong0.02spring1e3vLt.jpg" alt="prolong0.02spring1e3vLt" style="zoom:25%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.1spring1e3vLt.jpg" alt="0.1spring1e3vLt" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\expand0.1spring1e3vLt.jpg" alt="expand0.1spring1e3vLt" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\prolong0.1spring1e3vLt.jpg" alt="prolong0.1spring1e3vLt" style="zoom: 25%;" />

##### 球形方势阱

$$
V(x,y,z)=\left\{\begin{array}{lr}
k(\sqrt{x^2+y^2+z^2}-R),&r\geq R \\
0,&r<R
\end{array}\right.
$$



#### 分子链参数的影响

##### 键方差

| 模拟代号  | Parameter                                    |
| --------- | -------------------------------------------- |
| Normal    | `"bondLength": 1,"bondWiggleDistance": 0.05` |
| bondTight | `"bondLength": 1,"bondWiggleDistance": 0.01` |
| bondLoose | `"bondLength": 1,"bondWiggleDistance": 0.25` |

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.01spring1e3vLt.jpg" alt="0.01spring1e3vLt" style="zoom: 24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\bondParameter\figure\springBoneTight.jpg" alt="springBoneTight" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\bondParameter\figure\springBoneLoose.jpg" alt="springBoneLoose" style="zoom:24%;" />

可以看到，键长波动越大，最后的MSD波动越大。考虑到我们需要知道的力场参数和在MSD中可能和此项无关，可以考虑直接<u>**完全锁定键长**</u>。

##### 角度力

| 模拟代号   | Parameter  |
| ---------- | ---------- |
| Normal     | `"k": 1.5` |
| angleStiff | `"k": 7.5` |
| angleFlex  | `"k": 0.3` |

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.01spring1e3vLt.jpg" alt="0.01spring1e3vLt" style="zoom: 24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\bondParameter\figure\springAngleStiff.jpg" alt="springAngleStiff" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\bondParameter\figure\springAngleFlex.jpg" alt="springAngleFlex" style="zoom:24%;" />

观察MSD看不出太多区别，理论上Stiff会使得链变得更僵硬。

##### 交叠

| 模拟代号  | Parameter       |
| --------- | --------------- |
| Normal    | `"trunc": 1.0`  |
| crossLess | `"trunc": 10.0` |
| crossMore | `"trunc": 0.1`  |

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\0.01spring1e3vLt.jpg" alt="0.01spring1e3vLt" style="zoom: 24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\bondParameter\figure\springCrossMore.jpg" alt="springCrossMore" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\bondParameter\figure\springCrossLess.jpg" alt="springCrossLess" style="zoom:24%;" />

CrossMore代表的是，减少排斥力。可以观察到结出这个限制后，MSD波动显然增加变大。反之，CrossLess增加排斥力，使得波动变小。

#### 多次长时模拟

无力场情况下，考虑模拟次数增加对随机波动的影响：

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\100vLt No appled forcefield.jpg" alt="100vLt No appled forcefield" style="zoom: 24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\forceFieldVariableLangevin\figure\noForce1e3vLt.jpg" alt="noForce1e3vLt" style="zoom:24%;" />

<img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\longTime\figure\1000expNoForcefield10ps.jpg" alt="1000expNoForcefield10ps" style="zoom:24%;" /><img src="C:\Users\zhouquan\OneDrive\research\polychrom simulation\longTime\figure\1000expNoForcefield100ps.jpg" alt="1000expNoForcefield100ps" style="zoom:24%;" />
