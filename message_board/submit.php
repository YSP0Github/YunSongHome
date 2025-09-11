<?php

// 获取表单提交的数据
$name = $_POST['name'];
$email = $_POST['email'];
$message = $_POST['message'];

// 数据验证
if(empty($name) || empty($email) || empty($message)) {
  return ['status' => 'error', 'message' => '数据不完整'];
}

if(!filter_var($email, FILTER_VALIDATE_EMAIL)) {
  return ['status' => 'error', 'message' => '邮箱不合法'];  
}

// 保存到数据库
include 'db_connect.php'; // 包含数据库连接

$sql = "SELECT * FROM messages";

$result = mysqli_query($conn, $sql);

while($row = mysqli_fetch_assoc($result)) {
  echo $row['name'] . ": " . $row['email']. ": " . $row['message']; 
}

// 返回结果 
$result = [
  'status' => 'success',
  'message' => '留言成功'
];

echo json_encode($result);

?>