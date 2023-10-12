问题
- 当智能代理增加时， 使用主动声呐、雷达、和激光雷达来增强摄像头视觉的sota机器感知会遇到困难。而无处不在的热信号可以用来做补充
- 物体和环境无时无刻不散发着热信号，这会导致图片的'鬼影效应'

贡献
- 热辅助探测和测距解决了’鬼影问题‘，在黑暗条件下可以看到纹理和深度和其他物理属性
- 该工作带来了一种颠覆性的技术，这可以加速第四次工业革命

论文内容
- hue（色调）-> e， saturation（饱和度）-> T，brightness（亮度）-> X
- 哈达的工作波长由场景温度和大气透过率窗口决定，大气透过率窗口是可以在不被吸收或散射的情况下穿过大气层的波长范围
- $$S_{\alpha v} = e_{\alpha v}B_v(T_\alpha) + [1-e_{\alpha v}]X_{\alpha v}$$
	- 第一项是直接热辐射（无纹理），与发射率相关 ，==我暂时认为这是属于物体本身发出的辐射==
	- 携带纹理的第二项是==环境辐射经过物体表面散射后进入探测器==。反应了物体的光谱纹理对光谱辐射的影响，它考虑了影响反射或发射辐射的物体表面特性的变化
	- $_v$ 表示波数的依赖性
	- $B_v(T_\alpha)$ ：黑体在温度$T_\alpha$ 时的辐射， 受普朗克定律支配
 - 环境中其他物体$\beta$ 对 $\alpha$ 的热照度：$X_{\alpha v} = \sum_{\beta \neq \alpha}V_{\alpha \beta}S_{\beta v}$ 其中$V_{\alpha \beta}$ 是热照明系数。$\beta$ 是场景中大小有限的均匀紧凑物体
 - TeX degeneracy：$S_{\alpha v}$ 在T、e、X的联合变换下是不变的
 - 恢复纹理：打破TeX简并，并将光谱发射率$e_{\alpha v}$ 离散化到$e_v(m_\alpha)$ 中 ,材料库中包含了场景中所有可能的发射率
 - TeX-Net使用上面的公式来设计基于物理的损失，使用3D卷积神经网络来学习恢复X、T、e的空间光谱特征
### HADAR identifiability
- 物质差异由反射光谱之间的欧几里得距离决定（传统高光谱图片）
- 而哈达使用温度、纹理和发射率等多参数来识别物体。哈达可识别性定义为：可以从n个入射光子中检索到目标物质的最大香农信息，对应公式：$$I = log_2\Bigg\{1 + erf\bigg[\sqrt{\frac{Nd^2_0}{2(1+\gamma)}}\bigg]\Bigg\}$$ 
 
- $\gamma \equiv \gamma_1N + \gamma_0$ : 探测器的电子噪声功率经光子发射功率归一化后的值
- $d_0$ ：已知光谱发射率的两种材料之间的语义距离



### HADAR depth resolution
- 如果没有纹理信息，测距会变得不准确。而哈达可以把纹理回复到可以和光学图像相匹敌的程度，这样哈达测距和RGB立体视觉相当。


### Real-world HADAR preception
- 物理驱动感知而不是视觉驱动感知


## Methods
### TeX degeneracy
- 对于一个物体$\alpha$ ，它的光谱辐射$S_{\alpha v}$ 是不变的如果我们把它的物理属性$\{T_\alpha, e_{\alpha v}, X_{\alpha v}\}$ 改变成任意$\{T^{'}_{\alpha}, e^{'}_{\alpha v}, X^{'}_{\alpha v}\}$ 的话， 那么对应的：     ![[Pasted image 20230911194803.png]]   
	- v是波数， B是黑体辐射
- 具有不同的TeX属性但是具有相同的观测到的热信号$S_v$ 的物理状态被称为TeX简并

### TeX decomposition
- TeX-Net：同时使用空间信息和光谱信息
- 最小二乘估计和TeX-SGD（半局部分解）是非机器学习的baseline
### TeX vision and pseudo-TeX vision 
- 不同颜色代表不同类别，用HSV格式来表示TeX：H->e, S->T, V->X, 不同色调代表不同材料，饱和度对应温度，亮度对应纹理

### Eye safety restricts the scalability of LiDAR

### Monte Carlo path tracing
- HADAR数据库是长波红外立体高光谱数据库，由Blender Cycles渲染器中的普朗克定律和基尔霍夫定律合成。

HADAR estimation theroy
- 回答了需要多少光子来识别目标材料

### Prototype HADAR calibration and data collection 
- 通道数对应过滤器数

### Computational efficiency and deployablity
- TeX-Net 有50万个权重

### HADAR TeX vision algorithms
- TeX分解依赖于空间信息和光谱热特征，所以在网络中采用光谱和金字塔注意力层
- 由 $X_{\alpha v} = \sum_{\beta \neq \alpha}V_{\alpha \beta}S_{\beta \gamma}$ 必须被指定来保证逆映射的唯一性，所以要学习热照明系数V而不是纹理X， 因此它不是一个端到端的网络
	- $\alpha, \beta, \gamma$ 是物体的indices，v是波束
	- $S_{\beta \gamma}$ 被下采样到$S_{\alpha \gamma}$ 来近似k个最显著的环境物体
- 语料库和它的维度是网络的关键。
- 该网络可以用有监督训练也可以无监督训练



### 样机校准和数据收集
- prototype-1
	- detector : FLIRA325sc
	- 10个热红外过滤器检索光谱分辨率
	- 滤光轮上安装金镜监测探测器状态
	- Nicolet iS50 FTIR对滤光片的透过率进行表征
	- 标准黑体光源(EOII lnc DCN1000N7)标定相机的光谱响应曲线
	- 采集到的热立方体：**height \* width \* channel = 240 \* 320 \* 10** 通道数是过滤器的数量
- prototype-2
	- 带冷却碲化镉传感器的扫帚式高光谱成像仪： 提供256个光谱频段，价格超过100万美元


# Supplementary information
## Heat signal and heat cubes
- 恢复物体的**温度T**，**材料发射率 e** ， **纹理 X** 。
- 收集的是三维高光谱数据立方体$C_i(x, y, q)$ 
	- C : 探测器的输出数据
	- i ： 指的是探测器编号，即该热力方体由哪个探测器得到
	- q ：表示使用的滤波器
	- x,y : 探测器图像平面上的像素坐标
- 传感器收集小的物体元素辐射，而环境辐射只能通过物体元素的散射收集。![[Pasted image 20231009205151.png]]
- 目标物体$\alpha$ 的发射率为$e_\alpha$ ，到探测器的距离为z。环境中一小片区域$\beta$ 的 发射率为$e_\beta$ ,到目标物体的距离为$\rho$ ![[Pasted image 20231009210212.png]]
- 热照明细数V：在c中，探测器侧没有热辐射物体，所以v为0， d中，有一部分热辐射物体，所以$0 < V < 1$ ,e中，目标物体在探测器侧被环境包围，所以V为1![[Pasted image 20231009211309.png]]
- $\lambda_v$ : 输入到HADAR的热信号的平均光子数 

- 现有红外强度探测器主要基于热效应或光电二极管来响应光场强度
	- 前者吸收热辐射并将其转换成传感器的温度变化
	- 后者吸收热辐射并产生电子-空穴对

## 问题
- T,e,X具体如何转换到饱和度，亮度，色调？
- 温度和发射率为什么要分离？