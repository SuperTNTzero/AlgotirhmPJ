import matplotlib.pyplot as plt
import numpy as np

def reverse_complement(sequence):
    """生成DNA序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def hash_function(sequence, base_num=101, mod=10**9+7):
    """计算序列的哈希值
    使用简单的多项式哈希函数，适用于DNA序列
    优化版本：使用预计算的碱基映射和快速幂运算
    """
    # 预计算碱基映射，避免重复字典查找
    BASE_VAL = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    
    hash_val = 0
    for base in sequence:
        # 直接使用预计算的映射
        base_val = BASE_VAL.get(base, 0)
        hash_val = (hash_val * base_num + base_val) % mod
    return hash_val

def rolling_hash(prev_hash, prev_char, next_char, length, base_num=101, mod=10**9+7):
    """计算滚动哈希值，用于快速更新子序列的哈希值
    优化版本：使用预计算的碱基映射和幂值
    """
    # 预计算碱基映射，避免重复字典查找
    BASE_VAL = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    
    # 直接使用预计算的映射
    prev_val = BASE_VAL.get(prev_char, 0)
    next_val = BASE_VAL.get(next_char, 0)
    
    # 使用预计算的幂值
    highest_digit_val = (prev_val * pow(base_num, length-1, mod)) % mod
    
    # 移除最高位的贡献并添加新字符的贡献
    return ((prev_hash - highest_digit_val + mod) * base_num + next_val) % mod

def find_repeats_hash_dp(reference, query, min_length=3, max_length=None, use_parallel=True):
    """使用哈希和动态规划方法查找重复序列及其位置
    优化版本：添加提前终止条件和并行处理支持
    
    支持两种重复变异模式：
    1. 正向重复：某个序列s在reference的某个位置仅出现一次，但在query对应位置中连续出现了多次
       Sreference = −S1sS2−
       Squery = −S1sS3ss···ssS2−
    
    2. 反向重复：某个序列s在reference的某个位置仅出现一次，但其反向互补序列s'在query对应位置中连续出现了多次
       Sreference = −S1sS2−
       Squery = −S1sS3s′s′···s′s′S2−
    
    重复数量指的是重复片段在query中额外出现的重复次数（相比于reference额外出现的次数）
    
    参数:
        reference: 参考序列
        query: 查询序列
        min_length: 最小子序列长度
        max_length: 最大子序列长度，如果为None则使用两个序列的最小长度
        use_parallel: 是否使用并行处理
    """
    # 用于跟踪已经添加的序列，避免重复输出
    unique_sequences = set()
    import concurrent.futures
    import time
    
    results = []
    ref_len = len(reference)
    query_len = len(query)
    
    # 设置最大长度限制，避免处理过长的序列
    if max_length is None:
        max_length = min(ref_len, query_len)
    else:
        max_length = min(max_length, ref_len, query_len)
    
    # 如果序列太短，直接返回空结果
    if min_length > max_length:
        return results
    
    start_time = time.time()
    
    # 使用并行处理同时处理正向和反向序列
    if use_parallel and max_length > 20:  # 只有当序列足够长时才使用并行
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            # 提交正向序列处理任务
            forward_future = executor.submit(
                process_sequence, reference, query, False, [], min_length, max_length, unique_sequences
            )
            
            # 提交反向互补序列处理任务
            query_rev = reverse_complement(query)
            reverse_future = executor.submit(
                process_sequence, reference, query_rev, True, [], min_length, max_length, unique_sequences
            )
            
            # 获取结果
            forward_results = forward_future.result()
            reverse_results = reverse_future.result()
            
            # 合并结果
            results.extend(forward_results)
            results.extend(reverse_results)
    else:
        # 处理正向序列
        process_sequence(reference, query, False, results, min_length, max_length, unique_sequences)
        
        # 处理反向互补序列
        query_rev = reverse_complement(query)
        process_sequence(reference, query_rev, True, results, min_length, max_length, unique_sequences)
    
    # 按序列长度降序排序
    results.sort(key=lambda x: len(x['sequence']), reverse=True)
    
    end_time = time.time()
    print(f"查找重复序列耗时: {end_time - start_time:.2f} 秒")
    
    return results

def process_sequence(reference, query, is_reversed, results, min_length, max_length=None, unique_sequences=None):
    """处理序列查找重复子序列
    优化版本：使用滚动哈希和更高效的数据结构
    
    支持检测正向重复和反向重复变异：
    - 正向重复：某个序列s在reference的某个位置仅出现一次，但在query对应位置中连续出现了多次
    - 反向重复：某个序列s在reference的某个位置仅出现一次，但其反向互补序列s'在query对应位置中连续出现了多次
    
    返回找到的重复序列列表，支持并行处理
    
    参数:
        reference: 参考序列
        query: 查询序列
        is_reversed: 是否为反向重复
        results: 结果列表
        min_length: 最小子序列长度
        max_length: 最大子序列长度
        unique_sequences: 用于跟踪已处理序列的集合，避免重复输出
    """
    # 创建本地结果列表，支持并行处理
    local_results = [] if results is None else results
    
    ref_len = len(reference)
    query_len = len(query)
    base_num = 101
    mod = 10**9 + 7
    
    # 设置最大长度限制
    if max_length is None:
        max_length = min(ref_len, query_len)
    else:
        max_length = min(max_length, ref_len, query_len)
    
    # 预计算幂值，避免重复计算
    powers = [1]
    for i in range(1, max(ref_len, query_len)):
        powers.append((powers[-1] * base_num) % mod)
    
    # 为reference中所有可能的子序列创建哈希表
    ref_hash_map = {}
    
    # 对于每个可能的子序列长度
    for length in range(min_length, max_length + 1):
        # 提前终止条件：如果已经找到足够多的长序列，可以跳过短序列
        if len(local_results) > 100 and length < min_length + 5:
            break
            
        # 重置哈希表，只存储当前长度的子序列
        ref_hash_map.clear()
        
        # 计算reference中所有子序列的哈希值
        if ref_len >= length:
            curr_hash = hash_function(reference[:length], base_num, mod)
            
            if curr_hash not in ref_hash_map:
                ref_hash_map[curr_hash] = []
            ref_hash_map[curr_hash].append(0)
            
            # 使用滚动哈希计算其余子序列的哈希值
            for i in range(1, ref_len - length + 1):
                # 使用滚动哈希更新哈希值
                curr_hash = rolling_hash(curr_hash, reference[i-1], reference[i+length-1], length, base_num, mod)
                
                if curr_hash not in ref_hash_map:
                    ref_hash_map[curr_hash] = []
                ref_hash_map[curr_hash].append(i)
        
        # 计算query中所有子序列的哈希值并查找重复
        if query_len >= length:
            # 创建查询序列的哈希表，用于查找连续重复
            query_hash_map = {}
            
            # 计算query中第一个子序列的哈希值
            curr_hash = hash_function(query[:length], base_num, mod)
            
            # 将哈希值添加到查询哈希表
            if curr_hash not in query_hash_map:
                query_hash_map[curr_hash] = []
            query_hash_map[curr_hash].append(0)
            
            # 使用滚动哈希计算其余子序列的哈希值
            for i in range(1, query_len - length + 1):
                # 使用滚动哈希更新哈希值
                curr_hash = rolling_hash(curr_hash, query[i-1], query[i+length-1], length, base_num, mod)
                
                # 将哈希值添加到查询哈希表
                if curr_hash not in query_hash_map:
                    query_hash_map[curr_hash] = []
                query_hash_map[curr_hash].append(i)
            
            # 对于reference中的每个子序列，检查在query中是否有重复
            for hash_val, ref_positions in ref_hash_map.items():
                if hash_val in query_hash_map:
                    # 验证实际序列
                    for ref_pos in ref_positions:
                        ref_subseq = reference[ref_pos:ref_pos+length]
                        
                        # 如果序列已经处理过，则跳过
                        if unique_sequences is not None and ref_subseq in unique_sequences:
                            continue
                            
                        query_positions = []
                        
                        # 查找query中所有匹配的位置
                        for query_pos in query_hash_map[hash_val]:
                            query_subseq = query[query_pos:query_pos+length]
                            if ref_subseq == query_subseq:
                                query_positions.append(query_pos)
                        
                        # 如果找到重复（额外出现的次数大于0）
                        if len(query_positions) > 1:
                            # 计算额外重复次数（减去在reference中对应的一次出现）
                            repeat_count = len(query_positions) - 1
                            
                            # 添加重复结果
                            add_repeat_result(local_results, ref_subseq, [ref_pos], query_positions, is_reversed, repeat_count, length)
                            
                            # 将序列添加到已处理集合中
                            if unique_sequences is not None:
                                unique_sequences.add(ref_subseq)
    
    # 返回结果，支持并行处理
    return local_results

def add_result_if_not_contained(results, subseq, positions, is_reversed, length):
    """检查并添加结果，避免重复"""
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

def add_repeat_result(results, subseq, ref_positions, query_positions, is_reversed, repeat_count, length):
    """添加重复变异结果
    
    参数:
        results: 结果列表
        subseq: 重复序列
        ref_positions: 在reference中的位置列表
        query_positions: 在query中的位置列表
        is_reversed: 是否为反向重复
        repeat_count: 重复次数（额外出现的次数）
        length: 序列长度
    """
    # 检查是否已经有相同序列的结果
    for result in results:
        if result['sequence'] == subseq:
            # 序列已存在，不再添加
            return
    
    # 检查是否已经有更长的序列包含了这个子序列
    is_contained = False
    for result in results:
        if (result['reversed'] == is_reversed and 
            len(result['sequence']) > length and 
            all(p in result['ref_positions'] for p in ref_positions)):
            is_contained = True
            break
    
    if not is_contained:
        results.append({
            'sequence': subseq,
            'ref_positions': ref_positions,  # reference中的位置
            'query_positions': query_positions,  # query中的位置
            'repeat_count': repeat_count,  # 额外重复的次数
            'reversed': is_reversed
        })

def visualize_matches(reference, query, repeats, figsize=(12, 10), alpha=0.5, point_size=3, line_width=1.0):
    """可视化匹配结果，包括正向重复和反向重复变异
    
    参数:
        reference: 参考序列
        query: 查询序列
        repeats: 重复序列列表
        figsize: 图像尺寸，默认为(12, 10)
        alpha: 透明度，默认为0.5
        point_size: 点的大小，默认为3
        line_width: 线宽，默认为1.0
    """
    # 创建图形和子图
    fig, ax = plt.subplots(figsize=figsize)
    
    # 设置标题和轴标签
    ax.set_xlabel('Reference Sequence')
    ax.set_ylabel('Query Sequence')
    ax.set_title('DNA Sequence Repeat Variations')
    
    # 添加网格线
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # 绘制匹配点
    forward_x = []
    forward_y = []
    reverse_x = []
    reverse_y = []
    
    # 收集正向和反向匹配的坐标
    for repeat in repeats:
        seq_len = len(repeat['sequence'])
        
        # 检查是否有新格式的重复结果（包含ref_positions和query_positions）
        if 'ref_positions' in repeat and 'query_positions' in repeat:
            # 新格式的重复结果
            for ref_pos in repeat['ref_positions']:
                for query_pos in repeat['query_positions']:
                    if repeat['reversed']:
                        # 反向重复 - 绿色
                        if seq_len > 10:  # 增加长度阈值，减少线段数量
                            ax.plot([ref_pos, ref_pos + seq_len - 1], [query_pos, query_pos + seq_len - 1], 'g-', linewidth=line_width, alpha=alpha)
                        # 对于较短的序列，只添加端点，减少点的数量
                        if seq_len <= 10 or seq_len > 20:
                            reverse_x.append(ref_pos)
                            reverse_y.append(query_pos)
                            reverse_x.append(ref_pos + seq_len - 1)
                            reverse_y.append(query_pos + seq_len - 1)
                        else:
                            # 对于中等长度的序列，添加更多点以保持可见性
                            for i in range(0, seq_len, 2):  # 每隔2个碱基添加一个点
                                reverse_x.append(ref_pos + i)
                                reverse_y.append(query_pos + i)
                    else:
                        # 正向重复 - 红色
                        if seq_len > 10:  # 增加长度阈值，减少线段数量
                            ax.plot([ref_pos, ref_pos + seq_len - 1], [query_pos, query_pos + seq_len - 1], 'r-', linewidth=line_width, alpha=alpha)
                        # 对于较短的序列，只添加端点，减少点的数量
                        if seq_len <= 10 or seq_len > 20:
                            forward_x.append(ref_pos)
                            forward_y.append(query_pos)
                            forward_x.append(ref_pos + seq_len - 1)
                            forward_y.append(query_pos + seq_len - 1)
                        else:
                            # 对于中等长度的序列，添加更多点以保持可见性
                            for i in range(0, seq_len, 2):  # 每隔2个碱基添加一个点
                                forward_x.append(ref_pos + i)
                                forward_y.append(query_pos + i)
                    
                    # 在图上标注重复次数，但只为较长的序列添加标注，避免过度绘制
                    if 'repeat_count' in repeat and repeat['repeat_count'] > 0 and seq_len > 15:
                        ax.text(ref_pos + seq_len/2, query_pos - 5, f"重复{repeat['repeat_count']}次", 
                                fontsize=8, ha='center', va='bottom', 
                                color='green' if repeat['reversed'] else 'red')
        else:
            # 旧格式的重复结果（兼容性保留）
            for pos in repeat.get('positions', []):
                if repeat['reversed']:
                    # 反向互补匹配 - 绿色
                    # 对于反向互补，我们需要找到query中的对应位置
                    query_rev = reverse_complement(query)
                    subseq = repeat['sequence']
                    q_pos = query_rev.find(subseq)
                    if q_pos != -1:
                        # 只为较长的匹配添加线段，减少过度绘制
                        if seq_len > 10:  # 增加长度阈值，减少线段数量
                            # 添加线段
                            ax.plot([pos, pos + seq_len - 1], [q_pos, q_pos + seq_len - 1], 'g-', linewidth=line_width, alpha=alpha)
                        # 对于较短的序列，只添加端点，减少点的数量
                        if seq_len <= 10 or seq_len > 20:
                            reverse_x.append(pos)
                            reverse_y.append(q_pos)
                            reverse_x.append(pos + seq_len - 1)
                            reverse_y.append(q_pos + seq_len - 1)
                        else:
                            # 对于中等长度的序列，添加更多点以保持可见性
                            for i in range(0, seq_len, 2):  # 每隔2个碱基添加一个点
                                reverse_x.append(pos + i)
                                reverse_y.append(q_pos + i)
                else:
                    # 正向匹配 - 红色
                    # 找到query中的对应位置
                    subseq = repeat['sequence']
                    q_pos = query.find(subseq)
                    if q_pos != -1:
                        # 只为较长的匹配添加线段，减少过度绘制
                        if seq_len > 10:  # 增加长度阈值，减少线段数量
                            # 添加线段
                            ax.plot([pos, pos + seq_len - 1], [q_pos, q_pos + seq_len - 1], 'r-', linewidth=line_width, alpha=alpha)
                        # 对于较短的序列，只添加端点，减少点的数量
                        if seq_len <= 10 or seq_len > 20:
                            forward_x.append(pos)
                            forward_y.append(q_pos)
                            forward_x.append(pos + seq_len - 1)
                            forward_y.append(q_pos + seq_len - 1)
                        else:
                            # 对于中等长度的序列，添加更多点以保持可见性
                            for i in range(0, seq_len, 2):  # 每隔2个碱基添加一个点
                                forward_x.append(pos + i)
                                forward_y.append(q_pos + i)
    
    # 绘制正向匹配点（红色）
    forward_scatter = ax.scatter(forward_x, forward_y, color='red', s=point_size, alpha=alpha, label='Forward')
    # 绘制反向互补匹配点（绿色）
    reverse_scatter = ax.scatter(reverse_x, reverse_y, color='green', s=point_size, alpha=alpha, label='Reverse')
    
    # 设置单选按钮UI - 移动到左上角，避免与图形重叠
    from matplotlib.widgets import RadioButtons
    rax = plt.axes([0.02, 0.85, 0.15, 0.15], facecolor='lightgoldenrodyellow')
    radio = RadioButtons(rax, ('Both', 'Forward', 'Reverse'))
    
    # 创建图例标记 - 移动到左侧，避免与图形重叠
    legend_ax = plt.axes([0.02, 0.75, 0.15, 0.05], frameon=False)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    legend_ax.scatter([0.1], [0.7], color='red', s=30, label='Forward')
    legend_ax.scatter([0.1], [0.3], color='green', s=30, label='Reverse')
    legend_ax.text(0.2, 0.7, 'Forward', verticalalignment='center')
    legend_ax.text(0.2, 0.3, 'Reverse', verticalalignment='center')
    
    def update_visibility(label):
        # 更新线段和点的可见性
        if label == 'Both':
            for line in ax.lines:
                line.set_visible(True)
            for dot in ax.collections:
                dot.set_visible(True)
        elif label == 'Forward':
            # 显示正向匹配（红色）
            for i, line in enumerate(ax.lines):
                line.set_visible('r-' in str(line.get_color()))
            ax.collections[0].set_visible(True)
            ax.collections[1].set_visible(False)
        elif label == 'Reverse':
            # 显示反向互补匹配（绿色）
            for i, line in enumerate(ax.lines):
                line.set_visible('g-' in str(line.get_color()))
            ax.collections[0].set_visible(False)
            ax.collections[1].set_visible(True)
        fig.canvas.draw_idle()
    
    radio.on_clicked(update_visibility)
    
    # 显示图形
    plt.subplots_adjust(left=0.25, right=0.95, top=0.95, bottom=0.1)
    plt.show()

def main():
    # 获取输入序列
    file_path = input("请输入DNA序列文件路径（默认为'测试.txt'）：").strip()
    if not file_path:
        file_path = "测试.txt"
    
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            
            # 解析文件内容
            ref_start = content.find("ref:")
            query_start = content.find("query:")
            
            if ref_start == -1 or query_start == -1:
                print("文件格式错误，无法找到reference或query序列标记")
                return
            
            # 提取reference序列
            ref_content = content[ref_start + 4:query_start].strip()
            reference = ref_content.strip().upper()
            
            # 提取query序列
            query_content = content[query_start + 6:].strip()
            query = query_content.strip().upper()
            
            print(f"已从文件 {file_path} 读取序列：")
            print(f"Reference序列长度: {len(reference)}")
            print(f"Query序列长度: {len(query)}")
    except FileNotFoundError:
        print(f"文件 {file_path} 不存在")
        return
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return
    
    # 使用哈希和动态规划方法查找重复序列
    repeats = find_repeats_hash_dp(reference, query)
    
    # 输出结果
    print("\n找到的重复序列：")
    for i, repeat in enumerate(repeats[:10]):
        sequence = repeat['sequence']
        seq_length = len(sequence)
        is_reversed = repeat['reversed']
        
        # 检查是否有新格式的重复结果
        if 'ref_positions' in repeat and 'query_positions' in repeat:
            ref_pos = repeat['ref_positions'][0]  # 第一个位置
            repeat_count = repeat['repeat_count']  # 额外重复次数
            
            # 按照PPT要求输出信息
            print(f"重复 #{i+1}:")
            print(f"1. 第一次额外重复在reference的位置: {ref_pos}")
            print(f"2. 重复片段的长度: {seq_length}")
            print(f"3. 重复数量(额外出现的次数): {repeat_count}")
            print(f"4. 重复类型: {'反向重复' if is_reversed else '正向重复'}")
            print(f"重复序列: {sequence}")
        else:
            # 兼容旧格式
            positions = repeat.get('positions', [])
            count = repeat.get('count', 0) - 1  # 减去第一次出现
            
            if positions:
                print(f"重复 #{i+1}:")
                print(f"1. 第一次额外重复在reference的位置: {positions[0]}")
                print(f"2. 重复片段的长度: {seq_length}")
                print(f"3. 重复数量(额外出现的次数): {count}")
                print(f"4. 重复类型: {'反向重复' if is_reversed else '正向重复'}")
                print(f"重复序列: {sequence}")
        
        print()
    
    # 可视化匹配结果
    visualize_matches(reference, query, repeats)

if __name__ == "__main__":
    main()