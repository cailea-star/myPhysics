---
name: slides-html
description: 创建、编辑或审查以 slides.html 交付的静态学术与技术幻灯片。
---

# HTML 幻灯片

## 样式与结构

- 使用 4 空格缩进，遵守项目现有 HTML 与 CSS 格式。
- 必须优先复用已有变量、页面、内容块与布局。
- root.css 统一管理颜色、字号、高度与间距；优先调整现有变量。
- pages.css 负责页面结构及页眉页脚。
- blocks.css 负责内容块样式。
- layouts.css 负责行列排列、内容组合及块间距。
- 尺寸使用 --width-1pct 与 --height-1pct 规划。
- 不得使用行内样式或 <style> 覆盖公共样式。
- 每页使用一个 <section class="page ...">。

## 内容与资源

- 公式使用 MathJax v3。
- 图片使用相对路径，并保持原始比例。
- 公共 logo 存入 sildes/styles/。
- 页面专用图片存入该幻灯片目录的 assets/。
- 论文图片以 DOI 命名；同篇多图追加图号，并注明来源。
- 不得使用动画、过渡、Reveal.js 或主题系统。

## 工作流程

- 新增页面前先讨论内容与版型，取得批准后写入。
- 修改公共样式前明确影响范围，不得擅自扩大修改。
- 写入后检查溢出、对齐、公式渲染和资源路径。
- 导出使用 sildes/styles/export_pdf.py 输入.html 输出.pdf。
- 导出验证成功后才替换目标 PDF。
