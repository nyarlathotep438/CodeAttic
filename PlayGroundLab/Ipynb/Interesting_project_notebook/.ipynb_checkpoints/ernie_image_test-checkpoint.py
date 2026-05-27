# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 22:34:25 2026

@author: Admin
"""
import torch
from diffusers import ErnieImagePipeline

pipe = ErnieImagePipeline.from_pretrained(
    "Baidu/ERNIE-Image",
    torch_dtype=torch.bfloat16,
).to("cuda")

image = pipe(
    prompt="a cute cat sitting on a piano",
    height=1024,
    width=1024,
    num_inference_steps=20,
    guidance_scale=4.0,
    use_pe=True
).images[0]

image.save("test.png")
