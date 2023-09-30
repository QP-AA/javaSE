```c++
//泛化
template <class Key>
struct hash { };

//特化
template <>
struct  hash<char> {
	size_t operator() (char x) const {return x; }
};

template<>
struct hash<int> {
	size_t operator() (int x) const {return x; }
};
```
### partial specialization, 模板偏特化
```c++
//个数的偏
temlate<typename T, typename Alloc=.....>
class vector
{
...
};

template<typename Alloc=....>
class vector<bool, Alloc>
{

};

//范围的偏
template <typename T>
class C
{
...
};

template <typename T> //这里的T可以换成U V... 也就是说上面特化的T和下面的不是一个 
class C <T*>
{
...
};
```