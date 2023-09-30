```c++
template<typename T, template <typename T> class Container>
class XCLs {
private:
	Container<T> c;
public:
	...
};

template<typename T>
using Lst = list<T, allocator<T>>;

XCL<string, list> mylst1 //错误的
XCLs<string, Lst> mylst2 //正确的
```
由于list除过第一个参数之外还有一个分配器参数， 所以错误的是错误的， 正确的是正确的
### 关于c++标准库
啊？ 革命尚未成功， 同志仍需努力！

#### variadic template (since c++ 11)
数量不定的模板参数
```c++
void print() {
	//用于接收最后一次递归
}
template<typename T, typename... Types>
void print(const T& firstArg, const Types&... args) {
	cout << firstArg << endl;
	print(args...);
}
```
意思是把参数分为两部分，前面一个数量固定， 后面一包不固定， 数量任意。 注意...的位置

使用方法: `print(7.5, "hello", bitset<16>(377), 42);`
输出结果：
				7.5
				hello
				000000001011111
				43
使用`sizeof...(args)` 可以得到args的个数

#### auto (since c++ 11)
以前的版本：
```c++
lsit<string> c;
...
list<string>::iterator ite;
ite = find(c.begin(), c.end(), target);
```
使用auto后：
```c++
list<string> c;
...
auto ite = find(c.begin(), c.ed(), target);
```

#### ranged-base for (since c++ 11)
使用方法：
```c++
// decl 为变量 coll 为容器
for (decl : coll) {
	statement
}

//或者
for (int i : {2, 3, 4, 5, 6}) {
	cout << i << endl;
}

//再或者
vector<double> vec
...
for (auto elem : vec) {
	cout << elem << endl;
}
for (auto& elem : vec) {
	elem *= 3;
}
```