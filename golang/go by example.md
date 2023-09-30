go程序的执行：`go run xxx.go`
如果想将其编译成二进制文件：`go build xxx.go` 会生成一个.exe文件，输入`.\xxx.exe` 即可执行

对于变量来说：
	go 可以用:=的方式进行赋值，也可以用var来赋值。以下：
```go
	var a int = 3
	a := 3
```
两者都是把a的值赋成3

对于常量：
	用const来修饰。如`const a = "string"` . const可以出现在var出现的任何地方。

循环：
	`for`是go语言中唯一的循环形式
```go
	for i := 1; i <= 3; i++ {
		fmt.Println(i)
	}
```
if-else:
	条件选择语句，花括号必须。可以在if后面加一些声明语句，该声明可以在所有分支中使用。
```go
	if mun := 9; num < 0 {
		fmt.Println("negative")
	} else if num < 10 {
		fmt.Println("1-10")
	} else {
		fmt.Println(">10")
	}
```
	程序输出1-10

关于switch-case语句：
```go
		//一个基本的switch
		i := 3
		switch i {
			case 1:
			case 2:
			default:
		}
		//多个表达式
		 switch time.Now().Weekday() {
		    case time.Saturday, time.Sunday:
		        fmt.Println("It's the weekend")
		    default:
		        fmt.Println("It's a weekday")
		// 不带表达式
		 t := time.Now()
	    switch {
	    case t.Hour() < 12:
	        fmt.Println("It's before noon")
	    default:
	        fmt.Println("It's after noon")
	    }  
	    // 类型开关 i.(type)     
    }
```

空白标识符：_
```go
	m := make(map[string]int)
	_, pre := m["a"]
	// 由于m里面没有"a"，所以pre的值为false
```

对于range遍历，其可以用来迭代各种数据结构。在遍历数组，切片等数据结构时，第一个值为其下标。在遍历字符串时，返回的是字符的asc码

对于函数：
```go
	func plus(a,b int) int {
		return a + b
	}
```

入参的数据类型一致时，可以在最后一个入参后面加上对应的数据类型。后面的int是返回类型
如果为多返回值，则在返回类型内用括号将所有返回值类型都加上
```go
	func plusplus(a,b,c int) (int, int) {
		retun a, b
	}
```

变参函数：当传入的参数数量不确定时，可以如下方式定义函数
```go
	func sum(nums ...int) {
		fmt.Printlen(nums, " ")
	}
```

闭包：暂时不懂，看看例子：
```go
	package main

import "fmt"

func intSeq() func() int {
    i := 0
    return func() int {
        i++
        return i
    }
}

func main() {

    nextInt := intSeq()

    fmt.Println(nextInt())
    fmt.Println(nextInt())
    fmt.Println(nextInt())

    newInts := intSeq()
    fmt.Println(newInts())
}
```
函数执行结果为： 1， 2， 3,  1
闭包也可以是递归的，但是在定义闭包时，用类型化的var显示声明闭包：
```go
	var fib func(n int) int
	fib = func(n int) int {
		if n < 2 {
			return n
		}
		return fib(n - 1) + fib(n - 2)
	}
```
指针：
```go
	func a (i *int) {
		*i = 0  //对应值设为0
	}
	i := 1  // i = 1
	a(&i)  // &为取地址符， 此时 i值为0, &i为对应地址
```

字符串很难

结构体：
```go
type name struct {
	name  int
	age   int
}
```

方法：Go支持为结构体类型定义方法(methods)，在定义方法时，Go会自动处理值和指针之间的转换，如果想要避免在调用方法时产生一个拷贝，或者想让方法可以修改接受结构体的值，可以用指针来调用方法。