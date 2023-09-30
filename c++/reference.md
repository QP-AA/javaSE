```c++
int x = 0;
int* p = &x;  //个人理解 *p 指向x的地址 所以是&x 而不是 p = x
int& r = x; //有一个变量叫r 类型是reference to int， 直接输出r会输出x的值

```
引用不能再代表其他值， 如上代码 若`r = 2` 此时， x也会变成2。