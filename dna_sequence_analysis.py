import matplotlib.pyplot as plt
import numpy as np

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
    
    # 用于跟踪已经添加的序列，避免重复输出
    unique_sequences = set()
    
    # 正向搜索
    for i in range(query_len - min_length + 1):
        for length in range(min_length, query_len - i + 1):
            subseq = query[i:i+length]
            
            # 如果序列已经添加过，则跳过
            if subseq in unique_sequences:
                continue
                
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
                # 添加到已处理序列集合中
                unique_sequences.add(subseq)
    
    # 反向搜索
    query_rev = reverse_complement(query)
    for i in range(len(query_rev) - min_length + 1):
        for length in range(min_length, len(query_rev) - i + 1):
            subseq = query_rev[i:i+length]
            
            # 如果序列已经添加过，则跳过
            if subseq in unique_sequences:
                continue
                
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
                # 添加到已处理序列集合中
                unique_sequences.add(subseq)
    
    # 按序列长度降序排序
    results.sort(key=lambda x: len(x['sequence']), reverse=True)
    return results

def visualize_matches(reference, query, repeats, figsize=(12, 10), alpha=0.5, point_size=3, line_width=1.0):
    """可视化匹配结果
    
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
    ax.set_title('DNA Sequence Repeat Variations')
    ax.set_xlabel('Reference Sequence')
    ax.set_ylabel('Query Sequence')
    
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
        for pos in repeat['positions']:
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
    # 使用subplots_adjust代替tight_layout，避免与RadioButtons和legend_ax冲突
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
    
    # 可视化匹配结果
    visualize_matches(reference, query, repeats)

if __name__ == "__main__":
    main()