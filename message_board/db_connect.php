<?php
// db_connect.php

$db_host = "localhost";
$db_user = "root";
$db_pass = "password";
$db_name = "message_board";

$conn = mysqli_connect($db_host, $db_user, $db_pass, $db_name);

if (!$conn) {
  die("数据库连接失败: " . mysqli_connect_error());
}
?>