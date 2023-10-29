## 核心类和接口
- `DriverManager` 
	- `DriverManager.registerDriver(new Driver());` :会注册两次驱动,对应解决办法：`Class.forName("com.mysql.cj.jdbc.Driver");` ：通过反射直接调用
- `Connection`
	- `DriverManager.getConnection();` : ip+端口号， 具体库, username, password ![connectioin](jdbc02.png)
		- **url** : 语法：`jdbc:mysql://ip:Port/database` 
		- **Properties info** : 存储账号和密码， 需要创建对象
- `Statement`,`PreparedStatement`, `CallableStatement` 
	- `ResultSet` : 包含数据库返回的数据。内部包含游标，可以用`next` 进行移动，`next`返回值为true或者false。
	- `Statement` 存在的问题：
		- sql语句需要字符串拼接，比较麻烦
		- 只能拼接字符串类型，其他数据库类型无法处理
		- 可能发生注入攻击
	- `PreparedStatement` : 
		- 首先编写sql语句，不包含动态值部分的语句，动态值部分使用占位符 `?` 替代。 `？` 只能替代动态值
		- 创建PreparedStatement，并传入动态值
		- 动态值 占位符 赋值 `?`单独赋值
		- 发送sql并获取结果
- `Result`
	- `resultSet.getMetaData()` ：获取当前结果集中列的信息
## 使用路线
![jdbc](jdbc01.png)

- 具体操作步骤：
	- 以下为使用`statement`：
```java
// 注册驱动  
DriverManager.registerDriver(new Driver());  // 静态方法，可以通过类直接调用,cj.jdbc.Driver  由于注册两次驱动，所以采用下面的方法
Class.forName("com.mysql.cj.jdbc.Driver");   // 反射
  
// 获取连接  
Connection connection = DriverManager.getConnection("jdbc:mysql://127.0.0.1:3306/jdbc", "root", "12345");  
  
//创建statement  
Statement statement = connection.createStatement();  
  
// 发送sql语句， 获取返回结果  
String sql = "select * from t_user;";  
ResultSet resultSet = statement.executeQuery(sql);  
  
// 结果解析  
while (resultSet.next()) {  
    int id = resultSet.getInt("id");  
    String account = resultSet.getString("account");  
    String password = resultSet.getString("password");  
    String nickname = resultSet.getString("nickname");  
    System.out.println(id + "--" + account + "-- " + password + "--" + nickname);  
}  
  
// 关闭资源  
resultSet.close();  
statement.close();  
connection.close();
```

	- 以下为使用Preparedstatement进行查询
```java
// 其他步骤相同， 只是在创建和发送sql语句时略有不同
// ? 为动态值
String sql = "select * from t_user where account = ? and password = ?;";
// 创建PreparedStatement
PreparedStatement preparedStatement = connection.prepareStatement(sql);  
// 传入动态值， 1，2 表示动态值位置
preparedStatement.setObject(1, user_account);  
preparedStatement.setObject(2, user_password);  
  
// 发送并获取  
ResultSet resultSet = preparedStatement.executeQuery(); // 无需加参数sql
```
	- 进行更新：
```java
int ok = preparedStatement.executeUpdate(); // ok为影响的行数
```


## 连接池
