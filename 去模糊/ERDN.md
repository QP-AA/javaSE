![[Pasted image 20230814144210.png]]
![[Pasted image 20230814144227.png]]
$$F^a_{t+i}(p) = \Sigma_{k=1}^Kw_k*F_{t+i}(p+p_k+\Delta P_{t+i}(p)_k)$$
- $F_{t+i}$ 表示从相邻帧$I_{t+i}$ 所提取的特征。
- $w_k$ , $p_k$ 分别表示权重和预设的偏移量
- $\Delta P_{t+i}$ 可以从$F_t$ 和$F_{t+i}$ 的特征进行预测  
$$\Delta P^l_{t+i} = DB([F^l_{t+i}, F^l_t], (\Delta P^{l+1}_{t+i})^{\uparrow 2})$$ $$(F^a_{t+i})^l = f(DCN(F^l_{t+i}, \Delta P^l_{t+i}), ((F^a_{t+i})^{l+1})^{\uparrow 2})$$
- DCN 可变形卷积
- f 多个卷积层
- $\uparrow 2$ 代表两个上采样
- DB：扩张的空间金字塔



本文提出了一个观点：以前在每一个level使用很多卷积的方法会导致感受野的不一致，这会使得在较高层次预测的偏移高于较低层次的偏移。
本文提出的一个扩张的空间金字塔，连接几个具有不同膨胀速率的卷积层
$$C_i(F) = DConv_d(Conv_k(F))$$
- $DConv_d$ : 扩张卷积， 扩张速率为d， 卷积核大小为3x3， $Conv_k$ 代表传统卷积， 卷积核大小为k * k，![[Pasted image 20230815150255.png]]
