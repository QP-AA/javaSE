# 反射
- 创建对象的三种方式
	```java
// 类名.class  
Class clazz1 = Car.class;  
  
// 对象.class  
Class clazz2 = new Car().getClass();  
  
// class.forname  
Class clazz3 = Class.forName("com.jinwang.reflect.Car");  
  
  
// 实例化  无参构造
Car car = (Car)clazz3.getDeclaredConstructor().newInstance();

// 有参构造  
Constructor c1 = clazz.getConstructor(String.class, int.class, String.class);  
Car car = (Car) c1.newInstance("jieda", 10, "red");

// 获得属性
Field[] fields = clazz.getDeclaredFields(); 
```
