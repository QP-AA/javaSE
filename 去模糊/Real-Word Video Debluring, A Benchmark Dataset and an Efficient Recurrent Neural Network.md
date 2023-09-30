![[Pasted image 20230727102033.png]]![[Pasted image 20230727102053.png]]

- 从RDB-based RNN网络结构可以看出当前输入帧$I_t$ 将会被下采样，之后和最近的隐藏层$h_{t-1}$  进行一个concat操作得到$f_t^D$ , 接下来$f_t^D$ 经过一系列的RDB模块。得到$f_t^R$ , 而$h_t$ 则是在$f_t^R$ 的基础上再经过卷积->RDB->卷积

![[Pasted image 20230727103911.png]]


- GSA：用来提取和融合过去和未来帧的有效特征信息。GAP: global averaging pooling fusion, 以当前帧和相邻帧作为输入，然后过滤出有效分层特征，最后将所有特征concat并卷积得到输出。$$f^c_{t+i} = CAT(f_t, f_{t+i})$$ $$f^e_{t+i} = L(GAP(f_{t+i}^c))\bigotimes{P(f^c_{t+i})}$$ $$F_t = Conv(CAT(f_{t-2}^e, f_{t-1}^e, f^e_{t+1}, f^e_{t+2}, f_t))$$ 

- 在算法上的主要贡献：
	- 提出全局注意力机制GSA来选择包含更多时间信息的特征进行融合
	- 将残差密集模块引入到到RNN网络中进行空间信息的提取
- 主要优点：减少计算成本