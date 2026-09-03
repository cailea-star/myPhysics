---
inclusion: fileMatch
fileMatchPattern: '*.md'
---

##### `.md` 文件规范
1. 说明文档要尽可能简洁, 严格遵循下述格式
2.1 函数声明: `Eigen::VectorXd our_func(const Eigen::Ref<const Eigen::MatrixXd>& matF2D_A, const Eigen::Ref<const Eigen::VectorXd>& psiF1D_x)`
2.2 函数位置: `???-??? lines`
2.3 物理公式: $a_{ij} = \int \psi^*(r)A(r)\psi_j(r)dr$ (注意与代码逻辑对应)
2.3.1 实现步骤: 若算法复杂, 在物理公式下列出关键实现步骤
2.4 输入变量1: `A (Eigen::MatrixXd)`, 公式对应 $A$...
2.6 输出变量1: `a (double)`, 公式对应 $a_{ij}$ ...
2.7 函数数据流: `outer_func <- our_func <- inner_func`

```markdown
# 代码文件 ```file.cpp``` 说明

### `func_1`

##### 函数定义

- 函数位置: `???-??? lines`

```cpp
Eigen::VectorXd our_func(const Eigen::Ref<const Eigen::VectorXd>& xF1D_x);
```

##### 物理公式

- 物理含义: ...

$
y = x^2 ...
$

- 实现步骤 (若复杂):
  1. 步骤一: ...
  2. 步骤二: ...

##### 变量对应

- 输入变量1: `xF1D_x (Eigen::VectorXd)`, 公式对应 $x$
- 输出变量1: `yF1D (Eigen::VectorXd)`, 公式对应 $y$


##### 函数数据流
- 调用情况: `outer_func <- our_func <- inner_func`
```