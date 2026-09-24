---
name: slides-html
description: 使用项目现有 CSS 创建、编辑和审查静态 HTML 幻灯片。
---

# HTML 幻灯片

## 模板与样式

- 公共 CSS 默认只读；仅在用户明确授权修改 CSS 时改动。
- 必须复用现有页面、布局、组件及其 HTML 容器结构。
- 禁止新增专用 CSS、复制样式、行内样式或 `<style>` 覆盖。
- root.css 管参数，pages.css 管页面，layouts.css 管布局，blocks.css 管组件。
- 授权修改样式时，优先复用变量；尺寸使用 `--width-1pct`、`--height-1pct`。

## 内容与资源

- 使用 4 空格缩进；每页一个 `<section class="page ...">`。
- 公式使用 MathJax v3；保留批准的符号、指标和算符顺序。
- 图片使用相对路径、保持原始比例；页面资源存入本套 `assets/`。
- 公共 logo 复用 `sildes/styles/`；论文图片以 DOI 命名，多图追加图号，注明来源。
- 禁止动画、过渡、Reveal.js 和主题系统。

## 工作流程

- 新增页面或重排内容前，先给内容与版型草稿；批准后写入。
- 溢出优先通过已有布局、公式分行或拆页解决；不得擅自缩字号、删内容或改 CSS。
- 写入后检查实际渲染：公式、资源、对齐、溢出及元素重叠。
- 交付前核对修改范围，保留用户已有修改；未经授权不得改动公共 CSS。
- 导出使用 `sildes/styles/export_pdf.py 输入.html 输出.pdf`；验证成功后才替换目标 PDF。
