# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 01:27:38 2024

@author: 96328
"""
def find_min_number(remainders, divisors):
    """
    使用中国剩余定理（Chinese Remainder Theorem, CRT）的简化版本来解决“物不知数”问题。
    :param remainders: 余数列表，与divisors中的除数相对应。
    :param divisors: 除数列表，与remainders中的余数相对应。
    :return: 满足所有除数和余数条件的最小正整数。
    """
    # 计算所有除数的乘积（即模数的乘积）
    product_of_divisors = 1
    for divisor in divisors:
        product_of_divisors *= divisor

    # 初始化结果变量
    result = 0

    # 对于每一组除数和余数，计算其在CRT中的贡献
    for i, (remainder, divisor) in enumerate(zip(remainders, divisors)):
        # 计算当前除数对应的模数的逆元（在CRT中需要用到）
        # 这里我们使用扩展欧几里得算法来计算逆元
        # 但由于题目中的除数都是质数，我们可以直接使用费马小定理来计算逆元
        # 即 a^(p-1) ≡ 1 (mod p)，所以 a^(p-2) ≡ a^(-1) (mod p)
        # 对于质数p，a的逆元就是a^(p-2) % p
        inverse = pow(divisor, product_of_divisors // divisor - 2, product_of_divisors // divisor)

        # 计算CRT中的贡献并累加到结果中
        result += remainder * inverse * (product_of_divisors // divisor)

    # 由于可能存在负数结果，我们需要对结果取模以确保它是正数
    result %= product_of_divisors

    return result

# 定义题目中的除数和余数
divisors = [3, 5, 7]
remainders = [2, 3, 2]

# 调用函数并打印结果
answer = find_min_number(remainders, divisors)
print(f"满足所有条件的最小正整数是：{answer}")
