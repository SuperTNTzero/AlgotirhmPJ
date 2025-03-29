def reverse_complement(sequence):
    """生成DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def smith_waterman(reference, query, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
    """实现Smith-Waterman局部序列比对算法"""
    m, n = len(reference), len(query)
    # 初始化得分矩阵和回溯矩阵
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[None] * (n + 1) for _ in range(m + 1)]
    
    # 填充得分矩阵
    max_score = 0
    max_pos = (0, 0)
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # 计算匹配得分
            match = score_matrix[i-1][j-1] + (match_score if reference[i-1] == query[j-1] else mismatch_penalty)
            # 计算插入和删除得分
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            
            # 选择最大得分
            score_matrix[i][j] = max(0, match, delete, insert)
            
            # 更新回溯矩阵
            if score_matrix[i][j] == 0:
                traceback[i][j] = None
            elif score_matrix[i][j] == match:
                traceback[i][j] = 'match'
            elif score_matrix[i][j] == delete:
                traceback[i][j] = 'delete'
            else:
                traceback[i][j] = 'insert'
            
            # 更新最大得分位置
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)
    
    return max_score, max_pos, traceback

def find_repeats(reference, query, min_length=3):
    """查找重复序列及其位置"""
    results = []
    query_len = len(query)
    ref_len = len(reference)
    
    # 正向搜索
    for i in range(query_len - min_length + 1):
        for length in range(min_length, query_len - i + 1):
            subseq = query[i:i+length]
            positions = []
            pos = 0
            
            # 在reference中查找匹配
            while True:
                pos = reference.find(subseq, pos)
                if pos == -1:
                    break
                positions.append(pos)
                pos += 1
            
            if len(positions) > 1:  # 找到重复
                results.append({
                    'sequence': subseq,
                    'positions': positions,
                    'count': len(positions),
                    'reversed': False
                })
    
    # 反向搜索
    query_rev = reverse_complement(query)
    for i in range(len(query_rev) - min_length + 1):
        for length in range(min_length, len(query_rev) - i + 1):
            subseq = query_rev[i:i+length]
            positions = []
            pos = 0
            
            while True:
                pos = reference.find(subseq, pos)
                if pos == -1:
                    break
                positions.append(pos)
                pos += 1
            
            if len(positions) > 1:
                results.append({
                    'sequence': subseq,
                    'positions': positions,
                    'count': len(positions),
                    'reversed': True
                })
    
    # 按序列长度降序排序
    results.sort(key=lambda x: len(x['sequence']), reverse=True)
    return results

def main():
    # 获取输入序列
    print("请输入reference序列：")
    reference = input().strip().upper()
    print("请输入query序列：")
    query = input().strip().upper()
    
    # 使用Smith-Waterman算法进行序列比对
    score, pos, traceback = smith_waterman(reference, query)
    
    # 查找重复序列
    repeats = find_repeats(reference, query)
    
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