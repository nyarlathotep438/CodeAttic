# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:55:43 2026

@author: Admin
"""

import base64
from IPython.display import HTML
import mimetypes

def embed_file(path, width=None, height=None):
    mime, _ = mimetypes.guess_type(path)
    if mime is None:
        raise ValueError(f"无法识别文件类型: {path}")

    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()

    # 图片
    if mime.startswith("image/"):
        style = ""
        if width: style += f' width="{width}"'
        if height: style += f' height="{height}"'
        html = f'<img src="data:{mime};base64,{data}"{style}>'
        return HTML(html)

    # 音频
    if mime.startswith("audio/"):
        html = f'''
        <audio controls>
            <source src="data:{mime};base64,{data}">
        </audio>
        '''
        return HTML(html)

    # 视频
    if mime.startswith("video/"):
        html = f'''
        <video controls {f'width="{width}"' if width else ""}>
            <source src="data:{mime};base64,{data}">
        </video>
        '''
        return HTML(html)

    # 其他文件 → 提供下载链接
    html = f'<a download href="data:{mime};base64,{data}">下载 {path}</a>'
    return HTML(html)