def reverse_complement(sequence):
    """生成DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def hash_function(sequence, base_num=101):
    """计算序列的哈希值
    使用简单的多项式哈希函数，适用于DNA序列
    """
    hash_val = 0
    for base in sequence:
        # 将碱基映射为数字: A=1, C=2, G=3, T=4
        base_val = {'A': 1, 'C': 2, 'G': 3, 'T': 4}.get(base, 0)
        hash_val = (hash_val * base_num + base_val) % (10**9 + 7)
    return hash_val

def rolling_hash(prev_hash, prev_char, next_char, length, base_num=101, mod=10**9+7):
    """计算滚动哈希值，用于快速更新子序列的哈希值"""
    # 将碱基映射为数字
    prev_val = {'A': 1, 'C': 2, 'G': 3, 'T': 4}.get(prev_char, 0)
    next_val = {'A': 1, 'C': 2, 'G': 3, 'T': 4}.get(next_char, 0)
    
    # 移除最高位的贡献并添加新字符的贡献
    return ((prev_hash - prev_val * pow(base_num, length-1, mod)) * base_num + next_val) % mod

def find_repeats_hash_dp(reference, query, min_length=3):
    """使用哈希和动态规划方法查找重复序列及其位置"""
    results = []
    ref_len = len(reference)
    query_len = len(query)
    
    # 处理正向序列
    process_sequence(reference, query, False, results, min_length)
    
    # 处理反向互补序列
    query_rev = reverse_complement(query)
    process_sequence(reference, query_rev, True, results, min_length)
    
    # 按序列长度降序排序
    results.sort(key=lambda x: len(x['sequence']), reverse=True)
    return results

def process_sequence(reference, query, is_reversed, results, min_length):
    """处理序列查找重复子序列"""
    ref_len = len(reference)
    query_len = len(query)
    
    # 为reference中所有可能的子序列创建哈希表
    ref_hash_map = {}
    
    # 对于每个可能的子序列长度
    for length in range(min_length, min(ref_len, query_len) + 1):
        # 重置哈希表，只存储当前长度的子序列
        ref_hash_map.clear()
        
        # 计算reference中所有长度为length的子序列的哈希值
        for i in range(ref_len - length + 1):
            subseq = reference[i:i+length]
            h = hash_function(subseq)
            
            if h not in ref_hash_map:
                ref_hash_map[h] = []
            ref_hash_map[h].append(i)
        
        # 在query中查找匹配
        for i in range(query_len - length + 1):
            subseq = query[i:i+length]
            h = hash_function(subseq)
            
            # 如果哈希值匹配，验证实际序列
            if h in ref_hash_map and len(ref_hash_map[h]) > 1:
                positions = []
                for pos in ref_hash_map[h]:
                    if reference[pos:pos+length] == subseq:
                        positions.append(pos)
                
                if len(positions) > 1:  # 找到重复
                    # 检查是否已经有更长的序列包含了这个子序列
                    is_contained = False
                    for result in results:
                        if (result['reversed'] == is_reversed and 
                            len(result['sequence']) > length and 
                            all(p in result['positions'] for p in positions)):
                            is_contained = True
                            break
                    
                    if not is_contained:
                        results.append({
                            'sequence': subseq,
                            'positions': positions,
                            'count': len(positions),
                            'reversed': is_reversed
                        })

def main():
    # 获取输入序列
    print("请输入reference序列：")
    reference = input().strip().upper()
    print("请输入query序列：")
    query = input().strip().upper()
    
    # 使用哈希和动态规划方法查找重复序列
    repeats = find_repeats_hash_dp(reference, query)
    
    # 输出结果
    print("\n找到的重复序列：")
    for i, repeat in enumerate(repeats[:10]):
        print(f"序列: {repeat['sequence']} (长度: {len(repeat['sequence'])})")
        print(f"位置: {repeat['positions']}")
        print(f"重复次数: {repeat['count']}")
        print(f"是否逆转: {'是' if repeat['reversed'] else '否'}")
        print()

if __name__ == "__main__":
    main()