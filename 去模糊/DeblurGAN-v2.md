![[Pasted image 20230611171108.png]]

### fpn_inception
`self.enc0` : 先经过一个卷积层`nn.Conv2d(3,32,kernel_size=3, stride=2,   padding=padding, bias=False)` , 然后再批量归一化`self.bn = nn.BatchNorm2d(32, eps=0.001, momentum=0.1,  affine=True)` , 最后经过`nn.ReLU` 激活。
`self.enc1` : 
- 先经过一个卷积层`nn.Conv2d(32,32,kernel_size=3, stride=1,   padding=padding, bias=False)` , 然后再批量归一化`self.bn = nn.BatchNorm2d(32, eps=0.001, momentum=0.1,  affine=True)` , 然后经过`nn.ReLU` 激活。
- 先经过一个卷积层`nn.Conv2d(32,64,kernel_size=3, stride=1,   padding=1, bias=False)` , 然后再批量归一化`self.bn = nn.BatchNorm2d(32, eps=0.001, momentum=0.1,  affine=True)` , 然后经过`nn.ReLU` 激活。
- 最后经过一个最大池化层`nn.MaxPool2d(3,stride=2)` 
`self.enc2` : 与 `self.enc1` 类似，只是卷积层参数不同
总而言之，这些enc函数都是一些卷积+池化的操作。

