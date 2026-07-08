# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:46:30 2026

@author: Admin
"""

import base64
from IPython.display import HTML
import mimetypes

def embed_image(path, width=None):
    # 自动检测 MIME 类型（支持 webp）
    mime, _ = mimetypes.guess_type(path)
    if mime is None:
        raise ValueError(f"无法识别文件类型: {path}")

    # 读取文件并转 base64
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()

    # 构造 HTML
    style = f' width="{width}"' if width else ""
    html = f'<img src="data:{mime};base64,{data}"{style}>'

    return HTML(html)

# embed_image("figure.webp", width=400)
