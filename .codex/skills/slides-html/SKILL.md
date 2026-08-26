---
name: slides-html
description: "创建、编辑或审查以 slides.html 交付的静态学术与技术幻灯片；适用于固定 16:9 HTML/CSS、MathJax 和逐页审批，不用于 PowerPoint、Google Slides、Reveal.js 或动画演示。"
---

# HTML 幻灯片

## CSS 与结构

- 新建时必须复制 `assets/base.html` 和 `assets/slides.css`。
- 必须以 `assets/slides.css` 作为唯一格式基准。
- 不得使用行内样式或 `<style>` 覆盖模板格式。
- 新增样式必须先获批准，再同步回模板 CSS。
- 每页必须仅使用一个 `<section class="slide ...">`。

## 内容

- 公式必须使用 MathJax v3，不得使用公式截图。
- 行内公式使用 `\( ... \)`，独立公式使用 `\[ ... \]`。
- 必须定义必要符号，并保持整套记号一致。
- 图片必须使用相对路径，并保持原始纵横比。
- 不得使用动画、过渡、Reveal.js 或主题系统。

## 工作流程

1. 新增每页前，必须讨论目的、版型、文字、公式、图片和说明。
2. 必须逐页明确批准；大纲批准无效，每次仅可写入一页。
3. 实质改变内容或版型必须重新批准。
4. 写入后必须按模板画布渲染检查；实质性修复必须重新批准。
