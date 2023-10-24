## 核心类和接口
- `DriverManager` 
	- `DriverManager.registerDriver(new Driver());` :会注册两次驱动,对应解决办法：`Class.forName("com.mysql.cj.jdbc.Driver");` ：通过反射直接调用
- `Connection`
	- `DriverManager.getConnection();` : ip+端口号， 具体库, username, password ![connectioin](jdbc02.png)
		- **url** : 语法：`jdbc:mysql://ip:Port/database` 
		- **Properties info** : 存储账号和密码， 需要创建对象
- `Statement`,`PreparedStatement`, `CallableStatement` 
	- `ResultSet` : 包含数据库返回的数据。内部包含游标，可以用`next` 进行移动，`next`返回值为true或者false。
- `Result`
## 使用路线
![jdbc](jdbc01.png)

- 具体操作步骤：
```java
// 注册驱动  
DriverManager.registerDriver(new Driver());  // 静态方法，可以通过类直接调用,cj.jdbc.Driver  
  
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
